#!/usr/bin/env python
#
import os
import sys
from img_local_utilities import print_info, print_error, print_one_of_the_input_variables_not_set_error_and_exit, read_from_file_to_get_a_list


def parse_cog_tigr_cathfunfam_smart_supfam_input_gff_files(src_file, dest_file, submission_id, gold_id):
    input_gff_file, dataLoad_path = None, None
    gene_oid, protein_database_type, percent_identity, align_length, query_start, query_end, subj_start, subj_end, evalue, bit_score = None, None, None, None, None, None, None, None, None, None
    content_list = []

    if src_file and submission_id and gold_id:
        gbp_home = None

        if not os.path.exists(src_file):
            file_name = os.path.basename(src_file)
            src_file = "str(DM_ARCHIVE) + /" +  str(submission_id) + "/" + str(file_name)

        if os.path.exists(src_file):
            input_gff_file = src_file
        else:
            print_info('parse_cog_tigr_cathfunfam_smart_supfam_input_gff_files','unable to locate %s file from the original location or the restored file in the dm_path' %src_file)
            print_error('parse_cog_tigr_cathfunfam_smart_supfam_input_gff_files','a restore from jamo is required before proceeding with submission_id %s' %submission_id)

        content_list = read_from_file_to_get_a_list(input_gff_file)
        if content_list:
            if len(content_list) > 0:
                if os.path.exists(dest_file):
                    os.remove(dest_file)

                with open(dest_file,'w') as f:
                    for line in content_list:
                        subj_start, subj_end, evalue = None, None, None
                        line_list = line.split('\t')
                        #print "line_list %s" %line_list

                        attributes_field_list = line_list[-1].split(';')
                        #print "attributes_field_list %s" %attributes_field_list                        gene_oid = line_list[0]
                        gene_oid = line_list[0]
                        protein_database_type = line_list[2]
                        percent_identity = attributes_field_list[1].split('=')[1]
                        query_start = line_list[3]
                        query_end = line_list[4]
                        if dest_file.endswith("tigr.txt") or dest_file.endswith("pfam.txt"):
                            evalue = attributes_field_list[3].split('=')[1]
                            subj_start = attributes_field_list[4].split('=')[1]
                            subj_end = attributes_field_list[5].split('=')[1]
                        else:
                            subj_start = attributes_field_list[6].split('=')[1]
                            subj_end = attributes_field_list[7].split('=')[1]
                            evalue = attributes_field_list[4].split('=')[1]

                        bit_score = line_list[5]
                        align_length = attributes_field_list[2].split('=')[1]

                        new_line = None
                        if dest_file.endswith("cog.txt") or dest_file.endswith("pfam.txt"):
                            new_line = str(gene_oid)+'\t'+str(protein_database_type)+'\t'+str(percent_identity)+'\t'+str(align_length)+'\t'+str(query_start)+'\t'+str(query_end)+'\t'+str(subj_start)+'\t'+str(subj_end)+'\t'+str(evalue)+'\t'+str(bit_score)+'\n'
                            f.write(new_line)
                        elif dest_file.endswith("cathfunfam.txt") or dest_file.endswith("smart.txt") or dest_file.endswith("supfam.txt") or dest_file.endswith("tigr.txt"):
                            new_line = str(gene_oid)+'\t'+str(protein_database_type)+'\t'+str(percent_identity)+'\t'+str(query_start)+'\t'+str(query_end)+'\t'+str(subj_start)+'\t'+str(subj_end)+'\t'+str(evalue)+'\t'+str(bit_score)+'\t'+str(align_length)+'\n'
                            f.write(new_line)
            else:
                print_error('parse_cog_tigr_cathfunfam_smart_supfam_input_gff_files','something is wrong here, the gff file %s is empty' %input_gff_file)
    else:
        print_one_of_the_input_variables_not_set_error_and_exit('parse_cog_tigr_cathfunfam_smart_supfam_input_gff_files',src_file, submission_id, gold_id)

