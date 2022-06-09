from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio import Entrez
import argparse
import pathlib
from pathlib import Path
import os



def get_database_files(file_path):
    ''' Iterates over each directory, collecting all the database files
    '''
    list_of_folders = []
    list_of_files = []
    for path in file_path:
        directory_path = path.as_posix()
        list_of_folders.append(directory_path)
        for folder in list_of_folders:
            for root, dirs, all_files in os.walk(folder):
                for data_files in all_files:
                    files = os.path.join(root, data_files)
                    list_of_files.append(files)
    return list_of_files

def extract_from_database_files(files):
    ''' Extracts each scientific name and sequence from RIVM database files
    '''
    scientific_name = []
    scientific_name_list= []
    final_sequences = []
 
    counter = 0
    file_counter = 0
    complete_dict_for_kraken = {}
    for database_file in files:
        file_name = Path(database_file).stem
     
        file_counter += 1
        with open(database_file, 'rb') as data:
            always_print_def = False
            always_print_or = False
            contents = data.readlines()

            merged_contents_list = []
            for line in contents:
                line = line.decode("iso-8859-1")
                merged_contents_list.append(line)

            #convert contents to string to search for strings in whole file
            new_contents = ' '.join(merged_contents_list)
            if 'ORGANISM' in new_contents: 
                sequence_list_organism = []
                for line in contents:
                    line = line.decode("iso-8859-1")
                    if 'ORGANISM' in line:
                        scientific_name_organism = line.replace('ORGANISM', '').replace('    ', '')
                       
                        
                        #When line with organism contains Unkown instead of scientific name
                        if 'Unknown' in line:
                        
                            for line in contents:
                                line = line.decode("iso-8859-1")
                                if 'DEFINITION' in line:
                                    try:
                                        # Get genus and species name out of the line
                                        elements = line.split(' ')
                                        genus = elements[2]
                                        species = elements[3]
                                        scientific_name = genus + ' ' + species
                                    # Prevents error if line is empty
                                    except IndexError:
                                        scientific_name = 'No scientific name'
                                    # replaces unkown with scientific name 
                                    scientific_name_organism = scientific_name.replace(',', '')
                                    
                                    line.replace(scientific_name_organism, scientific_name)
                                    scientific_name_list.append(scientific_name_organism)
                              

                            
                    if "ORIGIN" in line:
                       always_print_or = True
                    if always_print_or:
                        origin = line.replace('ORIGIN', '').replace('//', '')
                        sequence = ''.join((x for x in origin if not x.isdigit()))
                        sequence = sequence.upper()
                        sequence = ''.join(sequence.strip('/n'))
                        sequence_organism = sequence.strip().replace(" ", "")
                        sequence_list_organism.append(sequence_organism)

                
                sequence_list = list(filter(None,sequence_list_organism))
                organism_sequence_list = "".join (str (element) for element in sequence_list)
               
                # combine scientific names and sequences in list
                filename_sequence_combo_list = []
             
                filename_sequence_combo_list.append(scientific_name_organism.strip())
                filename_sequence_combo_list.append(organism_sequence_list)
              
         
                
                # list with data to dictionary (filename=key and value=scientific name + sequence)
                complete_dict_for_kraken[file_name] = filename_sequence_combo_list

            

            elif 'DEFINITION' in new_contents:
                
                sequence_list_definition_def = []

                for line in contents:
                    line = line.decode("iso-8859-1")
                    if 'DEFINITION' in line:
                        if len(line) > 14:
                            elements = line.split(' ')
                          
                            genus = elements[2]
                            species = elements[3]
                            scientific_name_definition = genus + ' ' + species

                            scientific_name_list.append(scientific_name_definition)

                    if 'ORIGIN' in line:
                        always_print_def = True
                    if always_print_def:

                        #TODO haal de \r\n eruit .strip()
                        origin = line.replace('ORIGIN', '').replace('//', '')
                        sequence = ''.join((x for x in origin if not x.isdigit()))
                        sequence = sequence.upper()
                        sequence = ''.join(sequence.strip())
                        sequence_list_definition_def.append(sequence)

                sequence_list = list(filter(None,sequence_list_definition_def))
                sequence_definition = "".join (str (element) for element in sequence_list)

                # combine scientific names and sequences in list 
                filename_sequence_combo_list_def = []
            
                filename_sequence_combo_list_def.append(scientific_name_definition.strip())
                filename_sequence_combo_list_def.append(sequence_definition)
            
                # list with data to dictionary (filename=key and value=scientific name + sequence)
                complete_dict_for_kraken[file_name] = filename_sequence_combo_list_def
                final_sequences.append(sequence_definition)

    return complete_dict_for_kraken


def filter_dict(sequence_dict_incl_taxid):
    filtered_dictionary = {}
    for item in sequence_dict_incl_taxid:
        list_with_species_info2 = sequence_dict_incl_taxid[item]
        tax_id = list_with_species_info2[0]
        #filter if there is no taxid
        if len(tax_id) > 0:
            sequence_dict_incl_taxid[item]
            filetered_dictionary[item] = sequence_dict_incl_taxid[item]   
    return filtered_dictionary


def get_taxonomy_ID(complete_dict_for_kraken):
    ''' Search for the taxonomy ID based on the scientific name
    '''
    for item in complete_dict_for_kraken:
        list_with_species_info_incl_taxid = []
        list_of_species_information = complete_dict_for_kraken[item]
        species_name = list_of_species_information[0]
  
        Entrez.email = args.email
        handle = Entrez.esearch(db="Taxonomy", term=species_name)
        record = Entrez.read(handle)
        taxonomy_ID = record["IdList"]
    
        list_with_species_info_incl_taxid.append(taxonomy_ID)
        list_with_species_info_incl_taxid.append(list_of_species_information[0])
        list_with_species_info_incl_taxid.append(list_of_species_information[1])
        complete_dict_for_kraken[item] = list_with_species_info_incl_taxid
    return complete_dict_for_kraken

def create_multiFASTA(filtered_sequence_dictionary, output_dir):
    count = 0
    tax_id_string = '|kraken:taxid|' # creates header compatible with Kraken database
    counter_list = []
    for directory in output_dir:
        with open(str(directory)+"/combined_database.fasta", "w") as output_handle:
            for value in filtered_sequence_dictionary.values():
                #gets values from dictionary
                tax_id= value[0]
                taxonomy_ID = tax_id[0]
                scientific_name = value[1]
                count+=1 #counts to create unique sequence ID
                sequence = value[2]

                multifasta = '>sequence' + str(count) + tax_id_string + taxonomy_ID + "\t" + scientific_name + "\n" + sequence + "\n"
                output_handle.write(multifasta)
       

if __name__ == '__main__':
    argument_parser = argparse.ArgumentParser()
    argument_parser.add_argument('-db', '--database', type=pathlib.Path,
                        default=[], nargs='+', help='Path to folder with databases which need to be combined')
    argument_parser.add_argument('-e', '--email', type=str, help='Enter email from NCBI account')
    argument_parser.add_argument('-o', '--output', type=pathlib.Path,
                        default=[], nargs='+', help='Path to output directory')
    args = argument_parser.parse_args()
    extract_information = get_database_files(file_path=args.database)
    sequence_dictionary = extract_from_database_files(files=extract_information)
    sequence_dict_incl_taxid = get_taxonomy_ID(sequence_dictionary)
    filtered_seq_dict_incl_taxid = filter_dict(sequence_dict_incl_taxid)
    create_multiFASTA(filtered_seq_dict_incl_taxid, args.output)