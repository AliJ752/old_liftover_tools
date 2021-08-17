import pandas as pd
import numpy as np
import requests, sys
from tqdm import tqdm
import re

#filename = r"E:\Pipeline Shared\split_beds_for_analysis\20201014_PanHaemOnc_split_results\UCSC output\20201014_PanHaemOnc_01"
original_bed_filename = r"E:\Pipeline Shared\20201014_PanHaemOnc.bed"
UCSC_bed_filename = r"E:\Pipeline Shared\liftover_tests\20201014_PanHaemOnc_38_UCSC.bed"
cross_map_bed_filename = r"E:\Pipeline Shared\liftover_tests\20201014_PanHaemOnc_38_CrossMap.bed"
output_file = r"E:\Pipeline Shared\liftover_tests"

def read_in_bed(filepath):
    """
    Reads in a bed file directed to from filepath, 
    adds header columns

    :type filepath: string
    :param filepath: filepath to bed file
    """
    #table = pd.read_csv(filepath, sep='\t',header=None, 
    #skiprows=0,nrows=40) #rows limited for testing
    table = pd.read_csv(filepath, sep='\t',header=None) 
    table.rename(columns = {0: 'chr', 
                                 1: 'start',
                                    2: 'end', 
                                        3:'gene_symbol'}, inplace = True)
    return table

def calculate_interval_length(filepath):
    """
    Reads in a BED file and creates a 4th 
    column of the length of the interval.
    BED file is modified in place and returned.

    :type filepath: string
    :param filepath: filepath to bed file
    """
    bed_file = read_in_bed(filepath)
    bed_file['interval_lengths'] = (bed_file['end'] - bed_file['start'])

    return bed_file

def list_splits_and_unmapped(original_bed,liftover_bed):
    """
    Identifies those samples that have not 
    been mapped or have been split by crossmap.
    Takes the gene names from the bed files
    and sees if there are more following liftover 
    i.e split. If the gene name is no longer present 
    it has not been mapped. Mapped and unmapped 
    returned as lists.

    :type original_bed: dataframe
    :param original_bed: original bed with first four columns being 
                            chr,start,stop,gene_name

    :type liftover_bed: dataframe
    :param liftover_bed: liftover bed with first four columns being 
                            chr,start,stop,gene_name
    """
    original_gene_names = original_bed['gene_symbol'].tolist()
    lift_gene_names = liftover_bed['gene_symbol'].tolist()
    split = [gene for gene in original_gene_names if  
                    lift_gene_names.count(gene)>original_gene_names.count(gene)]
    unmapped = [gene for gene in original_gene_names if gene not in lift_gene_names]
    return unmapped,split
    
def compare_interval_size_nonsplit(original_bed,liftover_bed,unmapped_list,split_list):
    """
    Compares the interval sizes for those co-ordinates 
    that werent split by the liftover. Itterates through 
    the rest of the bed file via gene names, pulling out 
    interval values and comparing. Results not matching 
    returned in alt_interval_list.

    :type original_bed: dataframe
    :param original_bed: original bed with first five columns being 
                            chr,start,stop,gene_name,interval_lengths

    :type liftover_bed: dataframe
    :param liftover_bed: liftover bed with first four columns being 
                            chr,start,stop,gene_name,interval_lengths
    
    :type unmapped_list: list
    :param unmapped_list: list of unmapped co-ordinates returned from 
                            list_splits_and_unmapped

    :type split_list: list
    :param split_list: list of split co-ordinates returned from 
                            list_splits_and_unmapped

    """
    processed_original_bed = original_bed[~original_bed['gene_symbol'].isin(
        (unmapped_list + split_list))]
    processed_gene_names = processed_original_bed['gene_symbol'].tolist()
    alt_interval_sizes = []
    for gene in processed_gene_names:
        interval_in_original = processed_original_bed.loc[processed_original_bed[
            'gene_symbol'] == gene]['interval_lengths'].values
        interval_in_liftover  = liftover_bed.loc[liftover_bed[
            'gene_symbol'] == gene]['interval_lengths'].values
        if interval_in_original != interval_in_liftover:
            alt_interval_sizes.append(
                [gene,int(interval_in_original),int(interval_in_liftover)])
    return alt_interval_sizes

def compare_interval_size_split(original_bed,liftover_bed,unmapped_list,split_list):
    """
    Calculates and compares the interval size 
    of split co-ords.Identifies the single entry in
    the original bed and compares to the summed split
    entries.Returns altered interval list. 

    :type original_bed: dataframe
    :param original_bed: original bed with first five columns being 
                            chr,start,stop,gene_name,interval_lengths

    :type liftover_bed: dataframe
    :param liftover_bed: liftover bed with first four columns being 
                            chr,start,stop,gene_name,interval_lengths
    
    :type unmapped_list: list
    :param unmapped_list: list of unmapped co-ordinates returned from 
                            list_splits_and_unmapped

    :type split_list: list
    :param split_list: list of split co-ordinates returned from 
                            list_splits_and_unmapped    
    """
    altered_interval_sizes = []
    for gene in split_list:
        interval_in_original = original_bed.loc[
            original_bed['gene_symbol'] == gene]['interval_lengths'].values
        interval_in_liftover = liftover_bed[
            liftover_bed['gene_symbol'] == gene]['interval_lengths'].sum()
        if interval_in_original != interval_in_liftover:
            altered_interval_sizes.append([
                gene,int(interval_in_original),int(interval_in_liftover)])
    return altered_interval_sizes

def split_seperation(liftover_bed,split_list):
    """
    Calculates the difference between the end of 
    one split and start of the next. These are 
    returned as a nested list of gene and interval size.

    :type liftover_bed: dataframe
    :param liftover_bed: liftover bed with first four columns being 
                            chr,start,stop,gene_name,interval_lengths
    :type split_list: list
    :param split_list: list of split co-ordinates returned from 
                            list_splits_and_unmapped
    """
    seperation_values = []
    for gene in split_list:
        gene_table = liftover_bed[liftover_bed['gene_symbol'] == gene]
        end_of_first_split = int(gene_table.iloc[0:1, 2:3].values)
        start_of_next_split = int(gene_table.iloc[1:2, 1:2].values)
        seperation_values.append([gene,start_of_next_split-end_of_first_split])
    return seperation_values

def search_row_on_ensemble(bed_file,gene,server):
    """
    extracts interval chr,start and end data from a bed
    file for a single gene entry. Enters these into ensemble
    request. Returns search params to populate final table.

    :type original_bed: dataframe
    :param original_bed: original bed with first five columns being 
                            chr,start,stop,gene_name,interval_lengths
    :type gene: string
    :param gene: name of the gene to locate

    :type server: string
    :param server: url of reference build server
    """
    row = bed_file.loc[bed_file['gene_symbol'] == gene]
    #strips away to just int/char
    row_index = row.index.values[0]
    chr_number = (re.findall(r'[0-9]+|[X]|[Y]|[M]',str(row['chr'].values)))[0] 
    start = str(int(row['start'].values))
    end = str(int(row['end'].values))
    search = f"/overlap/region/human/{chr_number}:{start}-{end}?feature=gene;feature=transcript;feature=cds;feature=exon"
    search_params = [row_index,chr_number,start,end]
    #Make ensemble requests
    request_json = requests.get(server+search, headers={ "Content-Type" : "application/json"})
    return search_params,request_json

def check_features(original_bed,liftover_bed,unmapped_list,split_list,output_filename):
    """ 
    Itterates through the rows in both bed files using gene names. Searches 
    the ensemble API (37 and 38 as appropriate) using the chr,start,stop
    from each row in each bed file. Features and genes 
    are extracted from the returned json and compared. 
    Differences are compiled into a dataframe,
    returned, and saved as an excel file. 
 
 
    :type original_bed: dataframe
    :param original_bed: original bed with first five columns being 
                            chr,start,stop,gene_name,interval_lengths

    :type liftover_bed: dataframe
    :param liftover_bed: liftover bed with first four columns being 
                            chr,start,stop,gene_name,interval_lengths
    
    :type unmapped_list: list
    :param unmapped_list: list of unmapped co-ordinates returned from 
                            list_splits_and_unmapped

    :type split_list: list
    :param split_list: list of split co-ordinates returned from 
                            list_splits_and_unmapped       
    """
    #define respective rest API servers
    server_37,server_38 = "https://grch37.rest.ensembl.org","https://rest.ensembl.org"
    #get rid of all split and unmapped values 
    processed_original_bed = original_bed[~original_bed['gene_symbol'].isin((unmapped_list + split_list))]
    #define lists to capture all mismatches
    feature_mismatches = []
    gene_mismatches = []
    exon_mismatches = []
    print('Running feature comparison')
    for gene in tqdm(processed_original_bed['gene_symbol'].tolist()):
        # get row chr,start and stop values
        search_params37, build_37 = search_row_on_ensemble(processed_original_bed,gene,server_37)
        search_params38, build_38 = search_row_on_ensemble(liftover_bed,gene,server_38)
        returned_37_features_search,returned_38_features_search = [],[]
        returned_37_gene_search,returned_38_gene_search = [],[]
        returned_37_exon_search,returned_38_exon_search = [],[]
        #Itterate through returned json and seperate results into relevent lists
        for reference in [build_37,build_38]:
            if not reference.ok:
                reference.raise_for_status()
                sys.exit()
            decoded = reference.json() 
            for entry in decoded:
                if 'feature_type' in entry.keys():
                        try:
                            if reference == build_37:
                                if entry['feature_type'] not in returned_37_features_search:
                                    returned_37_features_search.append(entry['feature_type'])
                                if entry['feature_type'] == 'gene':
                                    returned_37_gene_search.append(entry['gene_id'])
                                if entry['feature_type'] == 'exon':
                                    returned_37_exon_search.append([entry['Parent'],entry['feature_type'],entry['rank']])
                                stored_37_entry_in_case_of_error = entry
                            else:
                                if entry['feature_type'] not in returned_38_features_search:
                                     returned_38_features_search.append(entry['feature_type'])
                                if entry['feature_type'] == 'gene':
                                    returned_38_gene_search.append(entry['gene_id'])     
                                if entry['feature_type'] == 'exon':
                                    returned_38_exon_search.append([entry['Parent'],entry['feature_type'],entry['rank']]) 
                            
                        #in case returned json has funny keys or is missing the keys
                        except (KeyError,TypeError) as e:
                            print(f'index {search_params37[0]} for gene {str(gene)} failed due to {e}')
                            print('Entry from build 37')
                            print(stored_37_entry_in_case_of_error)
                            print('Entry following liftover')
                            print(entry)
                            continue          
        #Adds to mismatch list if the two arent equal. 
        #Lists sorted so as lists with the same entries in different orders arent flagged as mismatched
        returned_37_features_search.sort()
        returned_38_features_search.sort()
        returned_37_gene_search.sort()
        returned_38_gene_search.sort() 
        returned_38_exon_search.sort()
        returned_38_exon_search.sort()

        if returned_37_gene_search != returned_38_gene_search:
            gene_mismatches.append([search_params37[0],gene,
            search_params37[1],search_params37[2],search_params37[3],
            search_params38[1],search_params38[2],search_params38[3],
            returned_37_gene_search,returned_38_gene_search,
            'NA','NA'])
        if returned_37_exon_search != returned_38_exon_search:
            missing_in_38 = [exon_details for exon_details in returned_37_exon_search if exon_details not in returned_38_exon_search]
            exon_mismatches.append([search_params37[0],gene,
            search_params37[1],search_params37[2],search_params37[3],
            search_params38[1],search_params38[2],search_params38[3],
            returned_37_gene_search,returned_38_gene_search,
            str(returned_37_exon_search),str(missing_in_38)])
        if (returned_37_features_search != returned_38_features_search):
            #the nested list thats added to the collated lists are just things I thought might be useful
            #chr_numbers are returned as lists for some reason which is why ive used an index.
            if not returned_37_features_search:
                continue
            missing_in_38 = [feature for feature in returned_37_features_search if feature not in returned_38_features_search]
            feature_mismatches.append([search_params37[0],gene,
            search_params37[1],search_params37[2],search_params37[3],
            search_params38[1],search_params38[2],search_params38[3],
            returned_37_gene_search,returned_38_gene_search,
            str(returned_37_features_search),str(missing_in_38)])
    #Makes a dataframe and adds column names for the data I chose to pull out.
    writer = pd.ExcelWriter(f'{output_file}{output_filename}',engine='xlsxwriter')
    dataframe_list = [gene_mismatches,exon_mismatches,feature_mismatches]
    for i in range(len(dataframe_list)):
        if i == 0:
            sheet_title = 'Gene mismatches'
        if i == 1:
            sheet_title = 'Exon mismatches'
        if i == 2:
            sheet_title = 'Feature mismatches'
        columns = ['row_index', 'specified_gene', 
                    'orig_chr_number','orig_start','orig_end',
                        'lift_chr_number','lift_start','lift_end',
                            'returned_37_genes','returned_38_genes',
                                'returned_37_features/exons','missing_features/exons_in_38']
        try:
            dataframe = dataframe_list[i]
            np.warnings.filterwarnings('ignore', category=np.VisibleDeprecationWarning)
            mismatch_df = pd.DataFrame(np.array(dataframe),
                                            columns=columns, dtype=object)
        except ValueError:
            print(f"No mismatches returned for {sheet_title} selection")
            mismatch_df = pd.DataFrame(columns=columns)
            mismatch_df = mismatch_df.fillna(0)

        mismatch_df.to_excel(writer,sheet_name=sheet_title,index=False)
        worksheet = writer.sheets[sheet_title]
        worksheet.autofilter('A1:M1')
    writer.save()

def main():
    original = calculate_interval_length(original_bed_filename)
    UCSC_bed = calculate_interval_length(UCSC_bed_filename)
    cross_map = calculate_interval_length(cross_map_bed_filename)
    cross_unmapped, cross_split = list_splits_and_unmapped(original,cross_map)
    #ucsc_unmapped, ucsc_split = list_splits_and_unmapped(original,UCSC_bed)
    #compare_interval_size_nonsplit(original,UCSC_bed,ucsc_unmapped,ucsc_split)
    #compare_interval_size_nonsplit(original,cross_map,cross_unmapped,cross_split)
    #compare_interval_size_split(original,cross_map,cross_unmapped,cross_split)
    #print(split_seperation(cross_map,cross_split))
    #check_features(original,UCSC_bed,ucsc_unmapped,ucsc_split,'ucsc_mismatch_summary.xlsx')
    check_features(original,cross_map,cross_unmapped,cross_split,'cross_mismatch_test.xlsx')
    print('Completed Successfully')


if __name__ == '__main__':
    main()