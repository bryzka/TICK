from typing import List

from Bio.Blast import NCBIWWW, NCBIXML
from Bio import SeqIO
from argparse import ArgumentParser

from config import *


def get_orgs(path: str) -> List[str]:
    """
    Gets lits of organisms to blast sequences against from a file
    :param path: path to file containing list of organisms
    :return: list of organisms
    """
    with open(path, 'r') as f:
        ret = f.readlines()
    return ret


if __name__ == '__main__':
    parser = ArgumentParser()
    parser.add_argument("input_file", help="Fasta file containing sequences that we want to check.")
    parser.add_argument("--config_file", '-c', '-config', default="config", help="Path to config file.")
    parser.add_argument("--email", help="email for blasting purposes")

    args = parser.parse_args()

    req, opt = config(args.config_file)

    if req["BLAST"] in {'WWW', "", "NCBI"}:
        blast = NCBIWWW
        blast.email = args.email
    else:
        raise NotImplementedError

    orgs = get_orgs(req["orgs"])
    fasta = SeqIO.read(args.input_file, "fasta")
    blast_ret = blast.qblast('blastn', 'nt', fasta, entrez_query=", ".join(orgs))
    blast_ret = NCBIXML.parse(blast_ret)
    for record in blast_ret:
        print("Record:")
        for alignment in record.alignments:
            for hsp in alignment.hsps:
                print("***Alignment***")
                print(f"Title: {alignment.title}")
                print(f"Seq: {hsp.sbjct}")
                print()
    print(blast_ret)


