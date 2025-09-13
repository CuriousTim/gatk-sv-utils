#define _POSIX_C_SOURCE 202405L // getline, strdup
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>

#include "htslib/khash.h"
#include "htslib/hts.h"
#include "htslib/vcf.h"

#define err(fmt, ...)                                           \
	do { fprintf(stderr, "error: " fmt "\n", __VA_ARGS__);  \
	exit(EXIT_FAILURE); }                                   \
	while(0)

KHASH_SET_INIT_STR(strset)
KHASH_MAP_INIT_STR(carriers, khash_t(strset)*)

khash_t(strset) *vidset;
khash_t(strset) *sidset;

void usage(FILE *fp)
{
	fprintf(fp, "usage: set_denovo_gt <invcf> <denovo> <outvcf>\n");
}

/**
 * Add string to hash set and return the string.
 */
const char *strset_add(khash_t(strset) *h, const char *s)
{
	const char *k = 0;
	khiter_t i = kh_get(strset, h, s);
	if (i == kh_end(h)) {
		k = strdup(s);
		if (k == 0)
			err("%s", "failed to duplicate string");
		int ret;
		kh_put(strset, h, k, &ret);
		if (ret == -1)
			err("%s", "failed to add string to hash set");
	} else {
		k = kh_key(h, i);
	}

	return k;
}

/**
 * Parse one line of the denovo variants file.
 */
void tokenize(char *line, const char **vid, const char **sid)
{
	char *delim = strchr(line, '\t');
	if (delim == 0)
		err("%s","denovo line has 0 fields");
	else if (line == delim)
		err("%s", "denovo VID is empty");


	if (*(delim + 1) == '\0')
		err("%s", "denovo SID is empty");

	*delim = '\0';
	*vid = strset_add(vidset, line);
	*sid = strset_add(sidset, delim + 1);
}

/**
 * Add one denovo variant to the hash map.
 */
void add_dnentry(khash_t(carriers) *h, const char *vid, const char *sid)
{
	int ret;
	khash_t(strset) *ss;
	khiter_t k = kh_put(carriers, h, vid, &ret);
	if (ret == -1) {
		err("%s", "failed to add VID to hash map");
	} else if (ret == 0) {
		ss = kh_val(h, k);
	} else {
		ss = kh_init(strset);
		if (ss == 0)
			err("%s", "failed to initialize hash set");
	}

	kh_put(strset, ss, sid, &ret);
	if (ret == -1)
		err("%s", "failed to add SID to hash set");

	kh_val(h, k) = ss;
}

void strset_free(khash_t(strset) *s)
{
	for (khiter_t i = kh_begin(s); i != kh_end(s); ++i) {
		if (kh_exist(s, i))
			free((char*)kh_key(s, i));
	}
	kh_destroy(strset, s);
}

void carriers_free(khash_t(carriers) *h)
{
	for (khiter_t i = kh_begin(h); i != kh_end(h); ++i) {
		if (kh_exist(h, i))
			kh_destroy(strset, kh_val(h, i));
	}
	kh_destroy(carriers, h);
}

/**
 * Load the denovos into a hash map.
 */
void load_denovos(const char *path, khash_t(carriers) *h)
{
	FILE *fp = fopen(path, "r");
	if (fp == 0)
		err("%s", "could not open denovos");

	char *line = 0;
	size_t linecap = 0;
	ssize_t linelen;
	const char *vid;
	const char *sid;
	while ((linelen = getline(&line, &linecap, fp)) > 0) {
		if (*(line + linelen - 1) == '\n')
			*(line + linelen - 1) = '\0';

		tokenize(line, &vid, &sid);
		add_dnentry(h, vid, sid);
	}
}

/**
 * Check that the input VCF header defines a String SVTYPE INFO tag.
 */
void check_for_info(bcf_hdr_t *hdr)
{
	int id = bcf_hdr_id2int(hdr, BCF_DT_ID, "SVTYPE");
	if (!bcf_hdr_idinfo_exists(hdr, BCF_HL_INFO, id)
			|| bcf_hdr_id2type(hdr, BCF_HL_INFO, id) != BCF_HT_STR) {
		err("%s\n", "SVTYPE is not in the INFO field");
	}
}

void subset_hdr_to_gt(bcf_hdr_t *hdr)
{
	int idx = bcf_hdr_id2int(hdr, BCF_DT_ID, "GT");
	if (idx == -1)
		err("%s", "GT is missing from header");
	for (int i = 0; i < hdr->n[BCF_DT_ID]; ++i) {
		if (i == idx)
			continue;
		if (bcf_hdr_idinfo_exists(hdr, BCF_HL_FMT, i))
			bcf_hdr_remove(hdr, BCF_HL_FMT, hdr->id[BCF_DT_ID][i].key);
	}
}

/**
 * Subset VCF record to GT FORMAT field.
 */
void subset_rec_to_gt(const bcf_hdr_t *hdr, bcf1_t *rec)
{
	int idx = bcf_hdr_id2int(hdr, BCF_DT_ID, "GT");
	for (int i = 0; i < hdr->n[BCF_DT_ID]; ++i) {
		if (i == idx)
			continue;
		if (bcf_hdr_idinfo_exists(hdr, BCF_HL_FMT, i))
			bcf_update_format(hdr, rec, hdr->id[BCF_DT_ID][i].key, 0, 0, BCF_HT_INT);
	}
}

/**
 * Update genotypes so all de novo carriers are alternate and all others are
 * reference.
 */
void update_genotypes(const bcf_hdr_t *hdr, bcf1_t *rec, khash_t(strset) *h)
{
	int nsample = bcf_hdr_nsamples(hdr);
	int32_t *gt_arr = 0;
	int ngt_arr = 0;
	int ngt = bcf_get_genotypes(hdr, rec, &gt_arr, &ngt_arr);
	int ploidy = ngt / nsample;
	if (ploidy != 2)
		err("%s", "sample ploidy is not 2");

	for (int i = 0; i < nsample; ++i) {
		int32_t *p = gt_arr + i * ploidy;
		bool is_carrier = kh_get(strset, h, hdr->id[BCF_DT_SAMPLE][i].key) != kh_end(h);
		if (p[0] == bcf_int32_vector_end || p[1] == bcf_int32_vector_end)
			err("%s", "samples do not all have the same ploidy");
		p[0] = bcf_gt_unphased(0);
		p[1] = bcf_gt_unphased(is_carrier ? 1 : 0);
	}
	bcf_update_genotypes(hdr, rec, gt_arr, ngt);
	free(gt_arr);
}

bool is_cnv(const bcf_hdr_t *hdr, bcf1_t *rec)
{
	char svtype[4] = {0, 0, 0, 0};
	int n = 4;
	int ret = bcf_get_info_values(hdr, rec, "SVTYPE", (void **)&svtype, &n, BCF_HT_STR);
	if (ret < 0) {
		printf("%d\n", ret);
		err("%s", "could not check SVTYPE of record");
	}

	bool tmp = strcmp(svtype, "CNV") == 0;

	return tmp;
}

int main(int argc, char *argv[])
{
	if (argc != 4) {
		usage(stderr);
		exit(2);
	}

	vidset = kh_init(strset);
	sidset = kh_init(strset);

	khash_t(carriers) *h = kh_init(carriers);
	load_denovos(argv[2], h);

	htsFile *infp = hts_open(argv[1], "r");
	if (!infp)
		err("%s", "failed to open input VCF");
	htsFile *outfp = hts_open(argv[3], "wz");
	if (!outfp)
		err("%s", "failed to open output VCF");

	bcf_hdr_t *hdr = bcf_hdr_read(infp);
	check_for_info(hdr);
	subset_hdr_to_gt(hdr);
	if (bcf_hdr_write(outfp, hdr) != 0)
		err("%s", "failed to write header to output VCF");
	bcf1_t *rec = bcf_init();
	while (bcf_read(infp, hdr, rec) == 0) {
		if (bcf_unpack(rec, BCF_UN_INFO))
			err("%s", "failed to parse VCF record");

		khiter_t k = kh_get(carriers, h, rec->d.id);
		if (k == kh_end(h))
			continue;

		if (is_cnv(hdr, rec))
			err("%s", "CNV de novo variants are not supported");
		subset_rec_to_gt(hdr, rec);
		update_genotypes(hdr, rec, kh_val(h, k));
		if (bcf_write(outfp, hdr, rec) != 0)
			err("%s", "failed to write record to output VCF");
	}
	bcf_destroy(rec);
	bcf_hdr_destroy(hdr);
	hts_close(infp);
	hts_close(outfp);
	carriers_free(h);
	strset_free(vidset);
	strset_free(sidset);
}
