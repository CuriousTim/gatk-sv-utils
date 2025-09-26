#define _POSIX_C_SOURCE 202405L // getline, strdup
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "htslib/khash.h"
#include "htslib/hts.h"
#include "htslib/vcf.h"

#define err(fmt, ...)                                           \
	do { fprintf(stderr, "error: " fmt "\n", __VA_ARGS__);  \
	exit(EXIT_FAILURE); }                                   \
	while(0)

struct vec {
	int *b;
	size_t n;
};

KHASH_MAP_INIT_STR(map, struct vec)

khash_t(map) *filters;

void vec_push(struct vec *v, int x)
{
	int *tmp = realloc(v->b, v->n + 1);
	if (!tmp)
		err("%s", "OOM");

	v->b = tmp;
	v->b[v->n++] = x;
}

void usage(FILE *fp)
{
	fprintf(fp, "usage: add_vcf_filters <invcf> <outvcf> <filters> [<new_headers>]\n");
}

void add_new_headers(bcf_hdr_t *hdr, const char *path)
{
	FILE *fp = fopen(path, "r");
	if (!fp)
		err("%s", "failed to open new headers file");
	char *line = 0;
	size_t linecap = 0;
	ssize_t linelen;
	int parse_len;
	while ((linelen = getline(&line, &linecap, fp)) > 0) {
		bcf_hrec_t *hrec = bcf_hdr_parse_line(hdr, line, &parse_len);
		if (!hrec) {
			if (parse_len == 0)
				err("header line does not start with '##': '%s'", line);
			else if (parse_len == -1)
				err("%s", "OOM");
			else
				err("malformed header: '%s'", line);
		}
		if (bcf_hdr_add_hrec(hdr, hrec) == -1)
			err("%s", "OOM");
	}

	free(line);
	fclose(fp);
}

void store_filter(const char *vid, char *s, bcf_hdr_t *hdr)
{
	khiter_t p = kh_get(map, filters, vid);
	struct vec f = {0, 0};
	if (p == kh_end(filters)) {
		char *tmp = strdup(vid);
		if (!tmp)
			err("%s", "OOM");
		int ret;
		p = kh_put(map, filters, tmp, &ret);
		if (ret == -1)
			err("%s", "OOM");
	} else {
		f = kh_val(filters, p);
	}

	char *i = s;
	char *j;
	while ((j = strchr(i, ';'))) {
		*j = '\0';
		int id = bcf_hdr_id2int(hdr, BCF_DT_ID, i);
		if (id == -1)
			err("'%s' is not defined in the header", i);
		vec_push(&f, id);
	}

	int id = bcf_hdr_id2int(hdr, BCF_DT_ID, i);
	if (id == -1)
		err("'%s' is not defined in the header", i);
	vec_push(&f, id);

	kh_val(filters, p) = f;
}

void read_new_filters(const char *path, bcf_hdr_t *hdr)
{
	FILE *fp = fopen(path, "r");
	if (!fp)
		err("%s", "failed to open filters file");

	char *line = 0;
	size_t linecap = 0;
	ssize_t linelen;
	while ((linelen = getline(&line, &linecap, fp)) > 0) {
		if (*(line + linelen - 1) == '\n')
			*(line + linelen - 1) = '\0';
		char *p = strchr(line, '\t');
		if (!p)
			continue;
		*p++ = '\0';
		if (strlen(line) == 0 || strlen(p) == 0)
			continue;

		store_filter(line, p, hdr);
	}

	free(line);
	fclose(fp);
}

void filters_free(void)
{
	for (khiter_t p = kh_begin(filters); p != kh_end(filters); ++p) {
		if (kh_exist(filters, p)) {
			free((void *)kh_key(filters, p));
			free(kh_val(filters, p).b);
		}
	}

	kh_destroy(map, filters);
}

int main(int argc, char *argv[])
{
	if (argc != 4 && argc != 5) {
		usage(stderr);
		exit(2);
	}

	htsFile *infp = hts_open(argv[1], "r");
	if (!infp)
		err("%s", "failed to open input VCF");
	bcf_hdr_t *hdr = bcf_hdr_read(infp);
	if (!hdr)
		err("%s", "failed to parse input VCF header");

	if (argc == 5)
		add_new_headers(hdr, argv[4]);

	filters = kh_init(map);
	if (!filters)
		err("%s", "OOM");
	read_new_filters(argv[3], hdr);

	htsFile *outfp = hts_open(argv[2], "wz");
	if (!outfp)
		err("%s", "failed to open output VCF");
	if (bcf_hdr_write(outfp, hdr) != 0)
		err("%s", "failed to write header to output VCF");

	bcf1_t *rec = bcf_init();
	while (bcf_read(infp, hdr, rec) == 0) {
		bcf_unpack(rec, BCF_UN_FLT);

		khiter_t p = kh_get(map, filters, rec->d.id);
		if (p != kh_end(filters)) {
			struct vec *v = &kh_val(filters, p);
			for (size_t i = 0; i < v->n; ++i) {
				bcf_add_filter(hdr, rec, v->b[i]);
			}
		}

		if (bcf_write(outfp, hdr, rec) == -1)
			err("%s", "failed to write VCF record");
	}
	bcf_destroy(rec);
	bcf_hdr_destroy(hdr);
	hts_close(infp);
	hts_close(outfp);
	filters_free();
}
