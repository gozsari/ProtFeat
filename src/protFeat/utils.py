import os, platform, shutil
import pathlib
path_to_folder = pathlib.Path(__file__).parent.resolve()
class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKCYAN = '\033[96m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'


def find_pssm_missing_proteins(fasta_dict, pssm_dir):
    set_missing_prots = set()
    set_prots_pssm_exists = set()
    for file in os.listdir(pssm_dir):
        protein_id = file.split(".")[0]
        set_prots_pssm_exists.add(protein_id)

    for protein_id in set_prots_pssm_exists:
        file = protein_id + ".pssm"
        flag = False
        sequence = ""
        with open(pssm_dir+"/"+file, "r") as fp:
            for line in fp:
                list_line = line.strip().split()
                if len(list_line) > 0:
                    if list_line[0] == '1':
                        flag = True
                if len(list_line) == 0:
                    flag = False
                if flag:
                    sequence += list_line[1]
        if protein_id in fasta_dict:
            if sequence != fasta_dict[protein_id]:
                set_missing_prots.add(protein_id)
    set_missing_prots = set_missing_prots.union(set(fasta_dict.keys()) - set_prots_pssm_exists)
    return list(set_missing_prots)


def read_fasta_to_dict(input_dir, fasta_file, place_protein_id):
    fasta_dict = dict()
    sequence = ""
    prot_id = ""
    with open("{}/{}.fasta".format(input_dir, fasta_file), "r") as fp:
        for line in fp:
            if line[0] == '>':
                if prot_id != "" and prot_id not in fasta_dict:
                    fasta_dict[prot_id] = sequence
                prot_id = line.strip().split("|")[place_protein_id]
                if place_protein_id == 0:
                    prot_id = prot_id[1:]
                if prot_id not in fasta_dict:
                    sequence = ""
            else:
                sequence += line.strip()
        fasta_dict[prot_id] = sequence
    fp.close()
    return fasta_dict

def form_single_fasta_files(list_proteins_no_pssm, fasta_dict):
    path_single_fastas = "{}/temp_folder/single_fastas".format(path_to_folder)

    if os.path.isdir(path_single_fastas) == False:
        os.mkdir(path_single_fastas)

    for prot_id in list_proteins_no_pssm:
        fw = open("{}/{}.fasta".format(path_single_fastas, prot_id), "w")
        fw.write(">sp|{}\n{}\n".format(prot_id, fasta_dict[prot_id]))
        fw.close()


def form_missing_pssm_files(pssm_dir):
    path_single_fastas = "{}/temp_folder/single_fastas".format(path_to_folder)
    path_blast = "{}/ncbi-blast".format(path_to_folder)

    for filename in os.listdir(path_single_fastas):
        prot_id = filename.split(".")[0].strip()
        single_fasta_prot = "{}/{}".format(path_single_fastas, filename)
        pssmfile = "{}/{}.pssm".format(pssm_dir, prot_id)

        if ("Linux" in str(platform.platform())):
            os.system("{}/psiblast -db {}/uniref50_db/uniref50.blastdb -evalue 0.001 -query {} "
                      "-out_ascii_pssm {}  -out {}/outfile -num_iterations 3 -comp_based_stats 1" \
                      .format(path_blast, path_blast, single_fasta_prot, pssmfile, path_blast))
        elif ("Darwin" in str(platform.platform())):
            os.system("{}/psiblastMAC -db {}/uniref50_db/uniref50.blastdb -evalue 0.001 -query {} "
                      "-out_ascii_pssm {}  -out {}/outfile -num_iterations 3 -comp_based_stats 1" \
                      .format(path_blast, path_blast, single_fasta_prot, pssmfile, path_blast))
    try:
        shutil.rmtree(path_single_fastas)
    except OSError as e:
        print("Error: %s - %s." % (e.filename, e.strerror))
def download_extract_ncbi_blast():
    import wget
    from zipfile import ZipFile
    remote_url = 'http://slpred.kansil.org/ncbi-blast.zip'
    local_file = '{}/ncbi-blast.zip'.format(path_to_folder)
    wget.download(remote_url, local_file)
    # Create a ZipFile Object and load sample.zip in it
    with ZipFile(local_file, 'r') as zipObj:
        # Extract all the contents of zip file in current directory
        zipObj.extractall('{}'.format(path_to_folder))
        zipObj.close()
    for ncbi_file in os.listdir('{}/ncbi-blast'.format(path_to_folder)):
        path_file_ncbi = '{}/ncbi-blast/{}'.format(path_to_folder, ncbi_file)
        os.chmod(path_file_ncbi, 0o777)

def copy_pssms(fasta_dict):
    import requests
    for prot_id in fasta_dict:
        remote_url = 'https://slpred.kansil.org/swissprot_pssms/{}.pssm'.format(prot_id)
        local_file = '{}/pssm_files/{}.pssm'.format(path_to_folder, prot_id)
        data = requests.get(remote_url)
        try:
            data.raise_for_status()
        except requests.exceptions.HTTPError:
            continue
        with open(local_file, 'wb') as file:
                file.write(data.content)

def copy_form_pssm_matrices(fasta_dict):
    pssm_dir = "{}/pssm_files".format(path_to_folder)
    list_proteins_no_pssm1 = find_pssm_missing_proteins(fasta_dict, pssm_dir)
    if len(list_proteins_no_pssm1) != 0:
        copy_pssms(list_proteins_no_pssm1)
    list_proteins_no_pssm2 = find_pssm_missing_proteins(fasta_dict, pssm_dir)
    if len(list_proteins_no_pssm2) != 0:
        form_single_fasta_files(list_proteins_no_pssm2, fasta_dict)
        form_missing_pssm_files(pssm_dir)


def edit_extracted_features_POSSUM(temp_output_file, output_file, fasta_dict):
    with open(temp_output_file, 'r') as fp:
        fw = open(output_file, 'w')
        for line, prot_id in zip(fp, fasta_dict):
            fw.write('{}'.format(prot_id))
            list_line = line.strip().split(',')
            for item in list_line:
                if item == 'nan' or item == 'inf' or item == '-inf':
                    fw.write('\t0')
                else:
                    fw.write('\t{}'.format(item))
            fw.write('\n')
        fw.close()
    fp.close()
    os.remove(temp_output_file)

def edit_extracted_features_iFeature(temp_output_file, output_file, place_protein_id):
    with open(temp_output_file, 'r') as fp:
        fw = open(output_file, 'w')
        for line in fp:
            if line[0] == '#':
                continue
            else:
                line_split = line.strip().split('\t')
                protein_id = line_split[0].split('|')[place_protein_id]
                fw.write(protein_id)
                for feature in line_split[1:]:
                    fw.write('\t{}'.format(feature))
                fw.write('\n')
        fw.close()
    fp.close()
    os.remove(temp_output_file)

def download_extract_iFeature_data():
    import wget
    from zipfile import ZipFile
    remote_url = 'https://slpred.kansil.org/iFeature/data.zip'
    local_file = '{}/iFeature/data.zip'.format(path_to_folder)
    wget.download(remote_url, local_file)
    # Create a ZipFile Object and load sample.zip in it
    with ZipFile(local_file, 'r') as zipObj:
        # Extract all the contents of zip file in current directory
        zipObj.extractall('{}/iFeature/'.format(path_to_folder))
        zipObj.close()
