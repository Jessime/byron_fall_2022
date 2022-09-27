# Computing GC Content
# https://rosalind.info/problems/gc/

# Given: At most 10 DNA strings in FASTA format (of length at most 1 kbp each).
# Return: The ID of the string having the highest GC-content, followed by the GC-content of that string.

# Sample
# >Rosalind_6404
# CCTGCGGAAGATCGGCACTAGAATAGCCAGAACCGTTTCTCTGAGGCTTCCGGCCTTCCC
# TCCCACTAATAATTCTGAGG
# >Rosalind_5959
# CCATCGGTAGCGCATCCTTAGTCCAATTAAGTCCCTATCCAGGCGCTCCGCCGAAGGTCT
# ATATCCATTTGTCAGCAGACACGC
# >Rosalind_0808
# CCACCCTCGTGGTATGGCTAGGCATTCAGGAACCGGAGAACGCTTCAGACCAGCCCGGAC
# TGGGAACCTGCGGGCAGTAGGTGGAAT

# Output
# Rosalind_0808
# 60.919540

import sys
from pathlib import Path


def split_before(iterable, pred, maxsplit=-1):
    """Yield lists of items from *iterable*, where each list ends just before
    an item for which callable *pred* returns ``True``:

    Note!! This function does NOT contain a bug
    """
    if maxsplit == 0:
        yield list(iterable)
        return

    buf = []
    it = iter(iterable)
    for item in it:
        if pred(item) and buf:
            yield buf
            if maxsplit == 1:
                yield [item] + list(it)
                return
            buf = []
            maxsplit -= 1
        buf.append(item)
    if buf:
        yield buf


def is_header(line: str):
    """Returns True if line in .fa file is header of new sequence."""
    return line.startswith(">")


def compute_highest_gc_content(fasta_data: list[str]):
    """For a fasta file, return the sequence name with the highest GC content,
    as well as that GC content, expressed as a percentage.
    """
    max_header = None
    max_gc = -1

    sequences = list(split_before(fasta_data, is_header))
    for seq_lines in sequences:
        header = seq_lines[0].strip(">")
        seq = "".join(seq_lines)
        gc_count = seq.count("G") + seq.count("C")
        current_gc = round(100 * gc_count / len(seq), 6)
        if current_gc > max_gc:
            max_gc = current_gc
            max_header = header

    return max_header, max_gc


def run(input_path: str):
    """Main entrypoint.

    Handles:
    1. Reading input
    2. Computing GC content
    3. Expressing output
    """
    fasta_data = Path(input_path).read_text().split()
    header, gc = compute_highest_gc_content(fasta_data)
    print(header)
    print(gc)


if __name__ == "__main__":
    run(sys.argv[1])
