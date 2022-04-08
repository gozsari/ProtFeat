from protFeat.feature_extracter import extract_protein_feature
import argparse
if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="arguments of feature extracter module")
    parser.add_argument('--pf', type=str, default="aac_pssm", help='protein feature in POSSUM or iFeature')
    parser.add_argument('--ppid', type=str, default=1, help='the place of protein id in fasta header')
    parser.add_argument('--inpf', type=str, default="input_folder", help='path to fasta file directory')
    parser.add_argument('--fname', type=str, default="sample", help='fasta file name')

    args = parser.parse_args()

    protein_feature = args.pf
    place_protein_id = int(args.ppid)
    input_folder = args.inpf
    fasta_file_name = args.fname


    extract_protein_feature(protein_feature, place_protein_id,
                            input_folder, fasta_file_name)