import argparse
from pathlib import Path
import re

CORRECT_NAME_RE = re.compile(
    r"(chr(X|Y|[1-9]|1[0-9]|2[0-2]))_(\d+)-(\d+)_(.+)_(DEL|DUP)_(__.+)\.jpg"
)
MISSED_NAME_RE = re.compile(
    r".+~~(.+_(DUP|DEL)_(chr(X|Y|[1-9]|1[0-9]|2[0-2]))_(\d+)_(\d+))~~(.+).jpg"
)


# The two path parsers return a tuple with
# (GD ID, chr, start, end, carrier)
def parse_correct_path(path):
    m = CORRECT_NAME_RE.fullmatch(path.name)
    if m is not None:
        return m.group(5, 1, 3, 4, 7)

    return m


def parse_missed_path(path):
    m = MISSED_NAME_RE.fullmatch(path.name)
    if m is not None:
        return m.group(1, 3, 5, 6, 7)

    return m


def parse_plots(path, parser, results=dict()):
    for plot_path in path.glob("*.jpg"):
        parts = parser(plot_path)
        if parts is not None:
            coords, carriers = results.get(parts[0], (parts[1:4], set()))
            carriers.add(parts[4])
            results[parts[0]] = (coords, carriers)

    return results


def parse_comparison_table(path, parsed_plots):
    results = dict()
    with open(path, mode="r", encoding="utf-8") as f:
        for line in f:
            fields = line.rstrip("\n").split("\t")
            carriers = results.get(fields[0], set())
            if len(fields[6]) > 0:
                carriers.update(fields[6].split(","))
            if len(fields[7]) > 0:
                carriers.update(fields[7].split(","))
            if len(fields[9]) > 0:
                carriers.update(fields[9].split(","))
            if len(fields[10]) > 0:
                _, tmp = parsed_plots.get(fields[5], (None, set()))
                carriers.update(tmp.intersection(fields[10].split(",")))
            results[fields[0]] = carriers

    return results


def write_new_cnvs(x, path):
    with open(path, mode="w", encoding="utf-8") as f:
        for k, v in x.items():
            coords = "\t".join(v[0])
            carriers = ",".join(v[1])
            f.write(f"{coords}\t{k}\t{carriers}\n")


def write_remove_calls(x, path):
    with open(path, mode="w", encoding="utf-8") as f:
        for k, v in x.items():
            for sid in v:
                f.write(f"{k}\t{sid}\n")


def parse_args():
    parser = argparse.ArgumentParser(description="Make GD revision tables")
    parser.add_argument(
        "correct",
        type=Path,
        metavar="<correct>",
        help="Directory of manually reviewed correct GD plots",
    )
    parser.add_argument(
        "missed",
        type=Path,
        metavar="<missed>",
        help="Directory of manually reviewed correct GD plots from VCF",
    )
    parser.add_argument(
        "comparisons",
        type=Path,
        metavar="<comparisons>",
        help="File with GD to VCF comparisons",
    )
    parser.add_argument(
        "new_cnvs",
        type=Path,
        metavar="<new-cnvs>",
        help="Where to write the new CNVs table",
    )
    parser.add_argument(
        "remove_calls",
        type=Path,
        metavar="<remove-calls>",
        help="Where to write the remove calls table",
    )

    return parser.parse_args()


def main():
    argv = parse_args()
    parsed_plots = parse_plots(argv.correct, parse_correct_path)
    parsed_plots = parse_plots(argv.missed, parse_missed_path, parsed_plots)
    parsed_compares = parse_comparison_table(argv.comparisons, parsed_plots)
    write_new_cnvs(parsed_plots, argv.new_cnvs)
    write_remove_calls(parsed_compares, argv.remove_calls)


if __name__ == "__main__":
    main()
