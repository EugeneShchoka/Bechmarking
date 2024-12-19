import gzip


def is_gzipped(file):
    return str(file).endswith(".gz") or str(file).endswith(".bgz")


def open_func(
        file,
        read_header=False,
        header_start="#",
        skip_rows=0,
        to_list=False,
        to_dict=False,
        sep="\t",
):
    """for row count: headers or lines with empty spaces are counted"""
    gzipped = is_gzipped(file)
    file_open = gzip.open if gzipped else open
    row_count = 0
    if to_dict:
        read_header = True
    header = None
    with file_open(file) as infile:
        for xline in infile:
            row_count += 1
            line = (
                xline.decode("utf-8").replace("\n", "")
                if gzipped
                else xline.replace("\n", "")
            )
            if (
                    line == ""
                    or (line.startswith(header_start) and read_header is False)
                    or row_count <= skip_rows
            ):
                continue
            if not to_dict:
                yield line if not to_list else line.split(sep)
            else:
                fields = line.replace(header_start, "").split(sep)
                if header is None:
                    header = fields
                    continue
                else:
                    yield dict(zip(header, fields))


def load_json(json_file, encoding="utf-8", parse_float=None):
    import json

    json_file = str(json_file)
    if is_gzipped(json_file):
        with gzip.open(json_file) as fh:
            return json.loads(
                fh.read().decode(encoding=encoding), parse_float=parse_float
            )
    else:
        with open(json_file) as fh:
            return json.load(fh, parse_float=parse_float)


def parse_json_lines(infile):
    import json as json

    for line in open_func(infile):
        yield json.loads(line)
