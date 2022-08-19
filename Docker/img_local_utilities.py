#!/usr/bin/env python
#

import sys
import os
import time
import shutil
import base64
import socket
import datetime
import re
import csv
import tarfile
import logging
import logging.handlers
import math
import datetime
import json
import hashlib
from collections import OrderedDict
#from bson import json_util
import subprocess
from subprocess import Popen, PIPE

# for some stupid reason, in an if block, i.e. 'if some_value:', python treats None and 0 ( 0.0 ) similarly and will
# not enter the if statement-- use expliit if some_value != None: to set for unset variables. Or it will not work for zero values.

def write_file_from_list(some_list,some_file,filter_string=None):
    # writes a list into a file-- one element per line.
    # the optional filter_string is used if it is desirable to
    # avoid writing lines w/matching string.

    if some_list and some_file:
        if os.path.exists(some_file):
            os.remove(some_file)

        with open(some_file, 'w') as f:
            for some_element in some_list:
                if filter_string:
                    if not filter_string in some_element :
                        line = str(some_element) + '\n'
                        f.write(str(line))
                else:
                    line = str(some_element) + '\n'
                    f.write(str(line))
    else:
        print_one_of_the_input_variables_not_set_error_and_exit(write_file_from_list,some_list,some_file)

def sort_list_into_odd_and_even_lists(some_list):
    odd_list, even_list = [], []

    if some_list:
        try:
            even_list = [x for x in some_list if x%2 ==0]
            odd_list = [x for x in some_list if x%2 !=0]

            if not even_list:
                even_list = []
            if not odd_list:
                odd_list = []

        except Exception as e:
            print_error('sort_list_into_odd_and_even_lists','sort_list_into_odd_and_even_lists error: %s' %e)
            sys.exit(1)

    else:
        print_input_variable_not_set_error_and_exit('sort_list_into_odd_and_even_lists','some_list')

    return odd_list, even_list

def normal_round(some_number):
    # round values  w/.5 or more to higher integer
    rounded_value = None

    if some_number != None:
        if some_number - math.floor(some_number) < 0.5:
            rounded_value = math.floor(some_number)
        else:
            rounded_value = math.ceil(some_number)
    else:
        print_input_variable_not_set_error_and_exit('normal_round','some_number')

    return rounded_value

def swap_key_value_of_dict(some_dictionary):
    new_dictionary = None

    if some_dictionary:
        if isinstance(some_dictionary, dict):
            if len(some_dictionary) > 0:
                new_dictionary =  dict([(value, key) for key, value in some_dictionary.items()])
            else:
                print_error('swap_key_value_of_dict','%s is empty' %some_dictionary)
        else:
            print_error('swap_key_value_of_dict','the invoking variable to this function is not a type dictionary')
    else:
        print_input_variable_not_set_error_and_exit('swap_key_value_of_dict','some_dictionary')

    return new_dictionary

def get_grouped_list_of_peer_sublists(some_list, delimiter):
    # each text base element (containing a common delimited character is partitioned up and sorted into another list
    # with each section (or sub-element from each element grouped together).  The need for this functionality
    # occurs when we can comparing a bunch of different lineages of gene_oid which becomes the input list
    # and we want to group all the same level of sublineages together.  Once we group the various sublineages together,
    # we can then count and determine the most common sublineages and predict the overall scaffold_oid lineage from
    # all the various gene_oid lineages.
    #
    # i.e.
    # 
    # list = [ Bacteria;Firmicutes;Bacilli;Bacillales;Bacillaceae;Siminovitchia;Siminovitchia farraginis,
    #          Bacteria;Firmicutes;Bacilli;Bacillales;Bacillaceae;Siminovitchia;Siminovitchia fordii,
    #          Bacteria;Firmicutes;Bacilli;Bacillales;Bacillaceae;Weizmannia;Weizmannia coagulans ]
    # 
    # grouped_list_of_peeer_sublists = [[Bacteria, Bacteria, Bacteria], [Firmicutes, Firmicutes, Firmicutes],
    #                                   [Bacilli, Bacilli, Bacilli], [Bacillaceae, Bacillaceae, Bacillaceae],
    #                                   [Siminovitchia, Siminovitchia, Weizmannia], 
    #                                   [Siminovitchia farraginis, Siminovitchia fordii, Weizmannia coagulans]]
    #
    # for a max of 8 'sublevels'
    grouped_list_of_peer_sublists = []

    if some_list and delimiter:
        anotherlist1, anotherlist2, anotherlist3, anotherlist4 = [], [], [], []
        anotherlist5, anotherlist6, anotherlist7, anotherlist8 = [], [], [], []

        for element in some_list:
            sublist = element.split(delimiter)
            no_sublist = len(sublist)

            for number in range(0, no_sublist):
                if number == 0:
                    anotherlist1.append(sublist[number])
                elif number == 1:
                    anotherlist2.append(sublist[number])
                elif number == 2:
                    anotherlist3.append(sublist[number])
                elif number == 3:
                    anotherlist4.append(sublist[number])
                elif number == 4:
                    anotherlist5.append(sublist[number])
                elif number == 5:
                    anotherlist6.append(sublist[number])
                elif number == 6:
                    anotherlist7.append(sublist[number])
                elif number == 7:
                    anotherlist8.append(sublist[number])

        grouped_list_of_peer_sublists = [anotherlist1, anotherlist2, anotherlist3, anotherlist4, \
                                          anotherlist5, anotherlist6, anotherlist7, anotherlist8]
        # remove empty sublists
        for list_i in grouped_list_of_peer_sublists:
            if len(list_i) == 0:
                grouped_list_of_peer_sublists.remove(list_i)
    else:
        print_one_of_the_input_variables_not_set_error_and_exit('get_grouped_list_of_peer_sublists',some_list, delimiter_char)

    return grouped_list_of_peer_sublists

def get_most_frequent_element_and_count_in_list(some_list):
    highest_count, count, no_of_copies = None, None, None
    most_frequent_sublist = []
    count_dict = {}

    if some_list:
        highest_count = 0
        for unique_element in set(some_list):
            count = list(some_list).count(unique_element)
            if count > highest_count:
                highest_count = count

            count_dict[ unique_element ] = count

        for element, count in count_dict.items():
            if count == highest_count:
                most_frequent_sublist.append(element)
    else:
        print_input_variable_not_set_error_and_exit('get_most_frequent_element_and_count_in_list','some_list')

    return list(set(most_frequent_sublist)), highest_count

def print_std_out_and_logfile(message_list,logging_object):
    # 2 types of message format accespted:
    #
    # (1) message = 'invoking function name: message body text'
    # (2) message = ''  (to print an empty line in both std out and the logfile)
    #

    if message_list and logging_object:
        for message in message_list:
            if ':' in message:
                message_list = message.split(':')

                if len(message_list) == 2:
                    print_info(message_list[0],message_list[1])
                    logging_object.info(message)

                elif len(message_list) > 2:
                    # there may be ':' in the mssage content-- therefore, lets extract the 
                    # function name and reconnect the message string.
                    invoking_function_name = message_list[0]
                    message_list.remove(message_list[0])
                    orig_message = ':'.join(message_list)

                    print_info(invoking_function_name, orig_message)
                    logging_object.info(message)

                else:
                    print("")
                    print_info('print_std_out_and_logfile','message format error')
                    print("")
                    print_info('print_std_out_and_logfile','the message is %s' %message)
                    print("")
                    print_info('print_std_out_and_logfile','message format:')
                    print_info('print_std_out_and_logfile','\t<invoking function name>:<message body')
                    logging_object.info("")
                    print("")
                    sys.exit(1)

            elif message == '':
                print("")
                logging_object.info("")

            else:
                print("")
                print_info('print_std_out_and_logfile','message format error')
                print_info('print_std_out_and_logfile','message format:')
                print_info('print_std_out_and_logfile','\t<invoking function name>:<message body')
                logging_object.info("")
                print("")
                sys.exit(1)

    else:
        print_one_of_the_input_variables_not_set_error_and_exit('print_std_out_and_logfile','message_list', 'logging_object')


def log_session(invoking_function_name):
    dataLoad_path, logfile = None, None

    if invoking_function_name:
        if os.environ.get('DATALOAD_HOME'):
            dataLoad_path = os.environ.get('DATALOAD_HOME')
            logfile = str(dataLoad_path) + '/logs/' + str(invoking_function_name).replace('.py','.log')
        else:
            print_not_set_error_and_exit('log_session','DATALOAD_HOME')

        log_format = ( '[%(asctime)s] %(message)s' )

        logging.basicConfig(level=logging.INFO,
                                  format=log_format,
                                  filename=(logfile))

        log = logging.getLogger(invoking_function_name)
        if log:
            log.info("")
            log.info("**********************************************************************************************")

            msg_list = [ '', 'log_session: new session', '']
            print_std_out_and_logfile(msg_list, log)

        else:
            print_not_set_error_and_exit('log_session','log')

    else:
        print_input_variable_not_set_error_and_exit('log_session','invoking_function_name')

    return log

def write_merfs_start_oid_file(accounting):
    dataLoad_path, merfs_start_oid_file, submission_id, taxon_oid = None, None, None, None

    if accounting:
        if accounting.get('dataLoad_path'):
            dataLoad_path = accounting[ 'dataLoad_path' ]
        else:
            print_not_set_error_and_exit('write_merfs_start_oid_file', 'dataLoad_path')
        if accounting.get('submission_id'):
            submission_id = accounting[ 'submission_id' ]
        else:
            print_not_set_error_and_exit('write_merfs_start_oid_file', 'submission_id')

        if accounting.get('taxon_oid'):
            taxon_oid = accounting[ 'taxon_oid'  ]
        else:
            print("")
            print_info('write_merfs_start_oid_file','no taxon_oid has been assigned to submission_id %s' %submission_id)
            print_info('write_merfs_start_oid_file','skipping the processing of this submission_id')
            return 'skip_processing'

        merfs_start_oid_file = str(dataLoad_path) + "/data/merfs/merfs_start_oid.txt"

        # before we do anything, check if the IMG/M ER processing flag is already set.
        hostname = socket.gethostname()
        if hostname:
            processing_flag = str(dataLoad_path) + '/IMG_M_ER_pipeline_current_running_on_' + str(hostname)
            if os.path.exists(processing_flag):
                with open(merfs_start_oid_file, 'w') as f:
                    f.write(str(taxon_oid))

                sys.exit(0)

            else:
                print_info('write_merfs_start_oid_file','since %s is missing, the IMG/M ER Pipeline' %processing_flag)
                print_info('write_merfs_start_oid_file','can only process submission_id(s) with this flag on the dataLoad dir')
                print_info('write_merfs_start_oid_file','or can/will corrupt pipeline data if submssions are processed without this file')
                print_path_does_not_exist_error_and_exit('write_merfs_start_oid_file',processing_flag)
        else:
            print_not_set_error_and_exit('write_merfs_start_oid_file', 'hostname')
    else:
        print_input_variable_not_set_error_and_exit, print_dict('write_merfs_start_oid_file', 'accounting')


def read_merfs_start_oid_file(merfs_start_oid_file):
    stored_taxon_oid = None
    lines = []

    if merfs_start_oid_file:
        if os.path.exists(merfs_start_oid_file):
            with open(merfs_start_oid_file, 'r') as f:
                lines = f.readlines()
                stored_taxon_oid = convert_some_value_into_more_usable_format(lines[0])

            if isinstance(stored_taxon_oid, int):
                return stored_taxon_oid

            else:
                print_error('read_merfs_start_oid_file','the file %s content %s is not consistent with the format of the merfs_start_oid_file' %(merfs_start_oid_file,lines))
        else:
            print_path_does_not_exist_error_and_exit('read_merfs_start_oid_file',merfs_start_oid_file)
    else:
        print_input_variable_not_set_error_and_exit('read_merfs_start_oid_file','merfs_start_oid_file')


def clean_null_None(d):
    if d:
        for k, v in list(d.items()):
            if d[k] == 'None':
                d[k] = ''
            elif d[k] == None:
                d[k] = ''

    else:
        print_input_variable_not_set_error_and_exit('clean_null_None','d')

    return d

def center_string(some_string):
    columns, center_string = None, None

    if some_string:
        try:
            # First determine width of terminal
            rows, columns = subprocess.check_output(['stty', 'size']).decode().split() 

        except:
            # Script running from cron
            columns = 80

        # format string
        center_string = some_string.center(int(columns))

    else:
        print_input_variable_not_set_error_and_exit('center_string','some_string')

    return center_string

def bold_string(some_string):
    if some_string:
        bold='\033[1m'
        unbold='\033[0m'

        bold_string = str(bold) + str(some_string) + str(unbold)
    else:
        print_input_variable_not_set_error_and_exit('bold_string','some_string')

    return bold_string

def recast_type_str_if_applicable(some_value):
    recasted_var = None

    if some_value != None:
        if some_value != 0:
            if isinstance(some_value , str):
                if some_value.isdigit():
                    try:
                        recasted_var = int(some_value)

                    except:
                        try:
                            recasted_var = float(some_value)
                        except:
                            recasted_var = some_value
                else:
                    recasted_var = some_value

            elif isinstance( some_value, int):
                recasted_var = int(some_value)

            elif isinstance( some_value, float):
                recasted_var = float(some_value)
        else:
            recasted_var = some_value
    else:
        print_input_variable_not_set_error_and_exit('recast_type_str_if_applicable','some_value')

    return recasted_var

def convert_some_value_into_more_usable_format(some_value):
    # convert type str -> int or float
    # bytes > str (if applicable) -> int or float
    # bytearray > str (if applicable) -> int or float
    some_value_new_type, some_intern_value = None, None

    if some_value != None:
        if isinstance(some_value, str):
            # the value could be '318' -> 318
            some_value_new_type = recast_type_str_if_applicable(some_value)

        elif isinstance(some_value, float):
            # the value could be 318.0 -> 318
            some_value_new_type = recast_type_str_if_applicable(some_value)

        elif isinstance(some_value_new_type, int):
            some_value_new_type = some_value
            
        elif isinstance(some_value, bytes) or isinstance(some_value, bytearray):
            some_intern_value = some_value.decode('utf-8')

            if isinstance(some_intern_value, str):
                some_value_new_type = recast_type_str_if_applicable(some_intern_value)

            elif isinstance(some_intern_value, int):
                some_value_new_type = int(some_intern_value)
            
            elif isinstance(some_intern_value, float):
                some_value_new_type = float(some_intern_value)

        elif isinstance(some_value, set) or isinstance(some_value, list):
            some_value_new_type = fix_type_of_set_or_list(some_value)

        else:
            some_value_new_type = some_value
    else:
        some_value_new_type = None

    return some_value_new_type

def fix_type_of_set_or_list(some_set_or_list):
    some_transformed_set_or_list = None

    if some_set_or_list:
        if isinstance(some_set_or_list, set):
            some_transformed_set_or_list = set([ convert_some_value_into_more_usable_format(x) for x in list(some_set_or_list) ])
        elif isinstance(some_set_or_list, list):
            some_transformed_set_or_list = [ convert_some_value_into_more_usable_format(x) for x in some_set_or_list ]
        else:
            some_transformed_set_or_list = some_set_or_list
    else:
        print_input_variable_not_set_error_and_exit('fix_type_of_set_or_list','some_set_or_list')

    return some_transformed_set_or_list

def join_list_into_a_comma_delimited_string(some_list):
    new_comma_delimitted_string, element_type = None, None
    new_string_list = []

    if some_list:
        for some_element in some_list:
            if isinstance(some_element, int):
                element_type = "type_int"

            new_string_list.append(str(some_element))

        if len(new_string_list) > 0:
            if element_type == "type_int":
                string_list = list(map(str, new_string_list))
                new_comma_delimitted_string = ','.join(string_list)
            else:
                new_comma_delimitted_string = ','.join(new_string_list)
    else:
        print_input_variable_not_set_error_and_exit('join_list_into_a_comma_delimited_string','some_list')

    return new_comma_delimitted_string

def get_md5sum_from_file(file_name):
    hash_md5, md5sum = None, None

    if file_name:
        if os.path.exists(file_name):
            print_info('get_md5sum_from_file','calculating the md5sum for file_name %s' %file_name)
            with open(file_name, "rb") as f:
                data = f.read()
                md5sum = hashlib.md5(data).hexdigest()
            if md5sum:
                print_info('get_md5sum_from_file','md5sum: %s' %md5sum)
            else:
                print_not_set_error_and_exit('get_md5sum_from_file','md5sum')
        else:
            print_path_does_not_exist_error_and_exit('get_md5sum_from_file',file_name)
    else:
        print_input_variable_not_set_error_and_exit('get_md5sum_from_file','file_name')
                
    return md5sum

def divide_and_conqueror_loading_into_sqlite3_file(conqueror_loading_function, taxon_oid, input_src_file_to_load, sqlite3_file, chunk_size = 100000, ignore_header = False):

    if sqlite3_file:
        if os.path.exists(sqlite3_file) and os.path.exists(input_src_file_to_load):
            cache_list = []
            line_number = 0

            for line in open(input_src_file_to_load):
                line_number += 1
                if line_number == 1:
                    cache_list.append(line) if not ignore_header else None
                else:
                    cache_list.append(line)
        
                if line_number%chunk_size == 0:
                    #conqueror_loading_function(cache_list)
                    conqueror_loading_function(cache_list, taxon_oid, sqlite3_file)
                    cache_list = []

            if cache_list:
                #conqueror_loading_function(cache_list)
                conqueror_loading_function(cache_list, taxon_oid, sqlite3_file)
        else:
            print_path_does_not_exist_error_and_exit('divide_and_conqueror_loading_into_sqlite3_file',sqlite3_file)
            print_path_does_not_exist_error_and_exit('divide_and_conqueror_loading_into_sqlite3_file',input_src_file_to_load)
    else:
        print_input_variable_not_set_error_and_exit('divide_and_conqueror_loading_into_sqlite3_file','sqlite3_file')

def sql_quote(some_string):
    quoted_string = None

    if some_string:
        if '`' in some_string:
            some_string = some_string.replace('`','\`')
        if '(' in some_string:
            some_string = some_string.replace('(','\(')
        if ')' in some_string:
            some_string = some_string.replace(')','\)')
        if ','  in some_string:
            some_string = some_string.replace(',','\,')
        if "'" in some_string:
            some_string = some_string.replace("'","\'")
        if '&' in some_string:
            some_string = some_string.replace('&','\&')
        if '[' in some_string:
            some_string = some_string.replace('[','\[')
        if ']' in some_string:
            some_string = some_string.replace(']','\]')
        if '-' in some_string:
            some_string = some_string.replace('-','\-')
        if ';' in some_string:
            some_string = some_string.replace(';','\;')
        if '~' in some_string:
            some_string = some_string.replace('~','\~')

        quoted_string = some_string

    else:
        print_input_variable_not_set_error_and_exit('sql_quote','some_string')

    return quoted_string

def write_config_file(some_dictionary, file_name):

    if some_dictionary and file_name:
        with open(file_name, 'w') as f:
            for key, value in sorted(some_dictionary.items()):    
                f.write('.%s %s\n'  %(key, value))

    else:
        print_one_of_the_input_variables_not_set_error_and_exit('write_config_file',some_dictionary,file_name)

def write_json_file(some_dictionary, json_file_name):

    if some_dictionary and json_file_name:
        with open(json_file_name, 'w') as f:
            f.write(json.dumps(some_dictionary, default=json_util.default))
    else:  
        print_one_of_the_input_variables_not_set_error_and_exit('write_json_file',some_dictionary,json_file_name)

def generate_taxon_links_for_img_merfs_data(accounting):
    dataLoad_path, submission_id, taxon_oid = None, None, None

    if accounting:
        if accounting.get('dataLoad_path'):
            dataLoad_path = accounting[ 'dataLoad_path' ]
        else:
            print_not_set_error_and_exit('generate_taxon_links_for_img_merfs_data','dataLoad_path')

        if accounting.get('submission_id'):
            submission_id = accounting[ 'submission_id' ]
        else:
            print_not_set_error_and_exit('generate_taxon_links_for_img_merfs_data','submission_id')

        if accounting.get('taxon_oid'):
            taxon_oid = accounting[ 'taxon_oid' ]
        else:
            print_not_set_error_and_exit('generate_taxon_links_for_img_merfs_data','taxon_oid')

        data_dir = str(dataLoad_path) + "/data/merfs/submissions"
        taxon_dir  = str(dataLoad_path) + "/data/merfs/taxons"

        for root, dirs, files in os.walk(taxon_dir):
            for file_name in files:
                os.remove(os.path.join(root, file_name))

        for root, dirs, files in os.walk(data_dir):
            for file_name in files:
                submission_file = str(data_dir) + "/" + str(file_name)
                taxon_file_name = os.path.basename(file_name).replace(str(submission_id), str(taxon_oid))
                taxon_file = str(taxon_dir) + "/" + str(taxon_file_name)
                os.symlink(submission_file, taxon_file)

    else:
        print_input_variable_not_set_error_and_exit('generate_taxon_links_for_img_merfs_data','accounting')

    return accounting


def print_info(the_invoking_functions_name, message_text, logging_object=None):
    # Standard print stmt listing:
    #    (1) calling function
    #    (2) text mssage
    #    (3) and optionally, write to the logging_object also

    if the_invoking_functions_name and message_text:
        print_string = "   " + str(the_invoking_functions_name) + ': INFO: ' + str(message_text)
        print(print_string)

        if logging_object:
            log_message = str(the_invoking_functions_name) + ': ' + str(message_text)
            logging_object.info(log_message)
    else:
        print_one_of_the_input_variables_not_set_error_and_exit('print_info',the_invoking_functions_name,message_text)

def print_error(the_invoking_functions_name, message_text, logging_object=None, phb_flag=None):
    # I added an optional flag of all the various print_error functions in img_local_utilities.py. 
    # When processing one task as part of the pipeline sequence of steps, it is desirable to exit
    # when some condition was not meet (i.e. print and exit). However if we are processing a large list
    # of taxon_oid (say 130,000) it would be desirable to print the issue rather than exit. 
    # If we exited on every issue on a large list, it would take forever to finish while simultaneously
    # generating negative attention from the pointy headed bosses "is it done, is it done?"
    # hense the acronym phb in phb_flag -phajek
    #
    # FYI: the feature was not added to any of the print_input_variable... functions below since
    #      before calling any functions, the input variable should be check to see if it is 
    #      populated before invoking any function. -phajek
    # 
    # Standard print ERROR stmt listing:
    #    (1) calling function
    #    (2) text mssage
    # action: 
    #    (1) print error messagw
    #    (2) exit(1)

    if the_invoking_functions_name and message_text:
        print("")
        print_string = " *  " + str(the_invoking_functions_name) + ': ERROR: ' + str(message_text)
        print(print_string)
        print("")

        if logging_object:
            log_message = str(the_invoking_functions_name) + ': ERROR: ' + message_text
            logging_object.info('')
            logging_object.info(log_message)
            logging_object.info('')

        if phb_flag == None:
            exit(1)
    else:
        print_one_of_the_input_variables_not_set_error_and_exit('print_error',the_invoking_functions_name,message_text)


def print_not_set_error_and_exit(the_invoking_functions_name, variable_name, logging_object=None, phb_flag=None):
    # I added an optional flag of all the various print_error functions in img_local_utilities.py. 
    # When processing one task as part of the pipeline sequence of steps, it is desirable to exit
    # when some condition was not meet (i.e. print and exit). However if we are processing a large list
    # of taxon_oid (say 130,000) it would be desirable to print the issue rather than exit. 
    # If we exited on every issue on a large list, it would take forever to finish while simultaneously
    # generating negative attention from the pointy headed bosses "is it done, is it done?"
    # hense the acronym phb in phb_flag -phajek
    #
    # Standard error message to print out when the variable is not set with 
    # any value.
    #
    # example print_not_set_error_and_exit('bold_string','some_string')

    if the_invoking_functions_name and variable_name:
        print_string1 = " * " + str(the_invoking_functions_name) + ": ERROR: the variable " + str(variable_name) +  " is not set with any value."
        print_string2 = " * " + str(the_invoking_functions_name) + ": ERROR: exiting."

        print("")
        print(print_string1)
        print(print_string2)
        print("")

        if phb_flag == None:
            sys.exit(1)

    else:
        print_one_of_the_input_variables_not_set_error_and_exit('print_not_set_error_and_exit',the_invoking_functions_name,variable_name)


def print_input_variable_not_set_error_and_exit(the_invoking_functions_name, variable_name):
    # a similar function as above except, it specified the unset variable to be the input variable
    # (variable used when invoking any given function).

    if the_invoking_functions_name and variable_name:
        print_string1 = " * " + str(the_invoking_functions_name) + ": ERROR: the input variable " + str(variable_name) +  " is not set with any value."
        print_string2 = " * " + str(the_invoking_functions_name) + ": ERROR: exiting."

        print("")
        print(print_string1)
        print(print_string2)
        print("")
        sys.exit(1)

    else:
        print_one_of_the_input_variables_not_set_error_and_exit('print_input_variable_not_set_error_and_exit',the_invoking_functions_name,variable_name)


def print_one_of_the_input_variables_not_set_error_and_exit(the_invoking_functions_name, *argv):

    if the_invoking_functions_name:
        print_string1 = " * " + str(the_invoking_functions_name) + ": ERROR:  one of the following input variables is not set:"

        print("")
        print(print_string1)
        count = 1
        for arg in argv:
            another_string = None

            if count == 1:
                another_string = " * " + str(the_invoking_functions_name) + ": ERROR:  the first input varable %s:" %arg

            elif count == 2:
                another_string = " * " + str(the_invoking_functions_name) + ": ERROR:  the seconds input varable %s:" %arg

            elif count == 3:
                another_string = " * " + str(the_invoking_functions_name) + ": ERROR:  the third input varable %s:" %arg

            elif count == 4:
                another_string = " * " + str(the_invoking_functions_name) + ": ERROR:  the fourth input varable %s:" %arg

            elif count == 5:
                another_string = " * " + str(the_invoking_functions_name) + ": ERROR:  the fifth input varable %s:" %arg

            elif count == 6:
                another_string = " * " + str(the_invoking_functions_name) + ": ERROR:  the sixth input varable %s:" %arg

            else:
                another_string = " * " + str(the_invoking_functions_name) + ": ERROR:  another input varable %s:" %arg

            print(another_string)
            count = count + 1

        print("")
        sys.exit(1)

    else:
        print_input_variable_not_set_error_and_exit('print_one_of_the_input_variables_not_set_error_and_exit','the_invoking_functions_name')


def print_one_of_the_variables_not_set_error_and_exit(the_invoking_functions_name, *argv):

    if the_invoking_functions_name:
        print_string1 = " * " + str(the_invoking_functions_name) + ": ERROR:  one of the following variables is not set:"

        print("")
        print(print_string1)
        count = 1
        for arg in argv:
            another_string = None

            if count == 1:
                another_string = " * " + str(the_invoking_functions_name) + ": ERROR:  the first varable %s:" %arg

            elif count == 2:
                another_string = " * " + str(the_invoking_functions_name) + ": ERROR:  the seconds varable %s:" %arg

            elif count == 3:
                another_string = " * " + str(the_invoking_functions_name) + ": ERROR:  the third varable %s:" %arg

            elif count == 4:
                another_string = " * " + str(the_invoking_functions_name) + ": ERROR:  the fourth varable %s:" %arg

            elif count == 5:
                another_string = " * " + str(the_invoking_functions_name) + ": ERROR:  the fifth varable %s:" %arg

            elif count == 6:
                another_string = " * " + str(the_invoking_functions_name) + ": ERROR:  the sixth varable %s:" %arg

            else:
                another_string = " * " + str(the_invoking_functions_name) + ": ERROR:  another varable %s:" %arg

            print(another_string)
            count = count + 1

        print("")
        sys.exit(1)

    else:
        print_input_variable_not_set_error_and_exit('print_one_of_the_input_variables_not_set_error_and_exit',the_invoking_functions_name)


def print_unable_to_connect_to_db(the_invoking_functions_name, database_name):

    if the_invoking_functions_name and database_name:
        print_string = " * " + str(the_invoking_functions_name) + ": ERROR: unable to connect to the %s database" %database_name
        print("")
        print(print_string)
        print("")
        sys.exit(1)

    else:
        print_one_of_the_input_variables_not_set_error_and_exit('print_unable_to_connect_to_db',the_invoking_functions_name,database_name)

def print_path_does_not_exist_error_and_exit(the_invoking_functions_name,file_name):

    if the_invoking_functions_name and file_name:
        print_string1 = " * " + str(the_invoking_functions_name) + ": ERROR: the file " + str(file_name) + " does not exists on the file system."

        if not os.path.exists(file_name):
            print("")
            print(print_string1)
            print("")
            sys.exit(1)
    else:
        print_one_of_the_input_variables_not_set_error_and_exit('print_path_does_not_exist_error_and_exit',the_invoking_functions_name,file_name)

def run_bash_command(cmd):
    stderr, tmp_output, cmd_output, formatted_output = None, None, None, None

    if cmd:

        #p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE).decode('utf-8'))
        p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, encoding='utf8')
        stdout, stderr = p.communicate()
        if stdout:
            cmd_output = stdout
        elif p.returncode != 0:
            print_info('run_bash_command', 'std out %s' %stdout)
            print_info('run_bash_command', 'std err %s' %stderr)
            print_info('run_bash_command', 'the return code for %s is %s' %(cmd,str(p.returncode)))

        #except subprocess.CalledProcessError as e:
        #except p.CalledProcessError as e:
        #    print("")
        #    print((bold_string("\t HALTING SCRIPTS EXECUTION:")))
        #    print((bold_string("\t\t ERROR FROM CMD: %s" %cmd)))
        #    print("")
        #    cmd_output = p.output.decode('utf-8')
        #    print(cmd_output)
        #    if isinstance(output,list):
        #        for some_error in output:
        #            error_list = str(some_error).split('\n')
        #            for error_line in error_list:
        #                print(("\t\t %s" %format(error_line)))
        #    else:
        #        error_list = str(output).split('\n')
        #        for error_line in error_list:
        #            print(("\t\t %s" %format(error_line)))
#
#            print("")
#            print("")
#            sys.exit(1)
        if cmd_output != None:
            #print("")
            #print("")
            #print("cmd_output %s" %cmd_output)
            if '\n' in cmd_output:
                output = cmd_output.split('\n')
                for one_line in output:
                    if one_line != None:
                        #print("output %s" %output)
                        ansi_list = cut_ansi_string_into_parts(one_line)
                        for ansi_item in ansi_list:
                            ansi_item = list(ansi_item)
                            if ansi_item[1] == 'Reset':
                                ansi_item[1] = None

                            if ansi_item[0]:
                                print(("\t\t  %s" %ansi_item[0]))
                            if ansi_item[1]:
                                print(("\t\t  %s" %ansi_item[1]))
                            if ansi_item[2]:
                                print(("\t\t  %s" %ansi_item[2]))
                            if ansi_item[3]:
                                print(("\t\t  %s" %ansi_item[3]))
            else:
                output = cmd_output
                ansi_list = cut_ansi_string_into_parts(output)

                for ansi_item in ansi_list:
                    ansi_item = list(ansi_item)
                    if ansi_item[1] == 'Reset':
                        ansi_item[1] = None

                    if ansi_item[0]:
                        print(("\t\t  %s" %ansi_item[0]))
                    if ansi_item[1]:
                        print(("\t\t  %s" %ansi_item[1]))
                    if ansi_item[2]:
                        print(("\t\t  %s" %ansi_item[2]))
                    if ansi_item[3]:
                        print(("\t\t  %s" %ansi_item[3]))
                    
            print("")
            print((bold_string("\t   end of child process output")))
            print("")
    else:
        print_input_variable_not_set_error_and_exit('run_bash_command','cmd')

    return cmd_output

def cut_ansi_string_into_parts(string_with_ansi_codes):
    #
    # Converts a string with embedded ANSI Color Codes and parses it to create a list of tuples describing pieces of the input string.
    #
    # :param string_with_ansi_codes:
    #
    # :return: [(sty, str, str, str), ...] A list of tuples. Each tuple has format: (text, text color, background color, effects)
    # 
    tuple_list, new_tuple_list = [], []

    if string_with_ansi_codes:
        color_codes_english = ['Black', 'Red', 'Green', 'Yellow', 'Blue', 'Magenta', 'Cyan', 'White', 'Green', 'Green', 'Reset']
        color_codes = ["30m", "31m", "32m", "33m", "34m", "35m", "36m", "37m", "94m", "92m", "0m"]
        effect_codes_english = ['Italic', 'Underline', 'Slow Blink', 'Rapid Blink', 'Crossed Out']
        effect_codes = ["3m", "4m", "5m", "6m", "9m"]
        background_codes = ["40m", "41m", "42m", "43m", "44m", "45m", "46m", "47m"]
        background_codes_english = ["Black", "Red", "Green", "Yellow", "Blue", "Magenta", "Cyan", "White"]

        ansi_codes = color_codes + effect_codes

        string_list = string_with_ansi_codes.split("\\u001b[")

        if (len(string_list)) == 1:
            string_list = string_with_ansi_codes.split("\033[")

        for teststring in string_list:
            if teststring == string_with_ansi_codes:
                tuple_list += [(teststring, None, None, None)]
                break
            if any(code in teststring for code in ansi_codes):
                static_string = None
                color_used = None
                effect_used = None
                background_used = None
                for color in range(0, len(color_codes)):
                    if teststring.startswith(color_codes[color]):
                        working_thread = teststring.split(color_codes[color])
                        ansi_strip = re.compile(r'\x1B[@-_][0-?]*[ -/]*[@-~]')
                        static_string = ansi_strip.sub('', working_thread[1])
                        color_used = color_codes_english[color]
                for effect in range(0, len(effect_codes)):
                    if teststring.startswith(effect_codes[effect]):
                        working_thread = teststring.split(effect_codes[effect])
                        ansi_strip = re.compile(r'\x1B[@-_][0-?]*[ -/]*[@-~]')
                        static_string = ansi_strip.sub('', working_thread[1])
                        effect_used = effect_codes_english[effect]
                for background in range(0, len(background_codes)):
                    if teststring.startswith(background_codes[background]):
                        working_thread = teststring.split(background_codes[background])
                        ansi_strip = re.compile(r'\x1B[@-_][0-?]*[ -/]*[@-~]')
                        static_string = ansi_strip.sub('', working_thread[1])
                        background_used = background_codes_english[background]
                try:
                    if not tuple_list[len(tuple_list) - 1][0]:
                        if not tuple_list[len(tuple_list) - 1][1] == None:
                            color_used = tuple_list[len(tuple_list) - 1][1]
                        if not tuple_list[len(tuple_list) - 1][2] == None:
                            background_used = tuple_list[len(tuple_list) - 1][2]
                        if not tuple_list[len(tuple_list) - 1][3] == None:
                            effect_used = tuple_list[len(tuple_list) - 1][3]
                        tuple_list += [(static_string, color_used, background_used, effect_used)]
                    else:
                        tuple_list += [(static_string, color_used, background_used, effect_used)]
                except Exception:
                    tuple_list += [(static_string, color_used, background_used, effect_used)]


        for x in range(0, len(tuple_list)):
            if tuple_list[x][0]:
                new_tuple_list += [(tuple_list[x][0], tuple_list[x][1], tuple_list[x][2], tuple_list[x][3])]

    #else:
    #    print_input_variable_not_set_error_and_exit('cut_ansi_string_into_parts', 'string_with_ansi_codes')

    return new_tuple_list



def print_dict(some_dictionary,dictionary_name):
    # example: to print out the contents of the dictionary accounting
    # print_dict(accounting, 'accounting')i
    another_dict = {}

    if some_dictionary and dictionary_name:
        if some_dictionary == {}:

            return

        elif len(some_dictionary) > 0:
            #print("print_dict function: type(some_dictionary) %s" %type(some_dictionary))
            print("")
            title_string = bold_string("              * Content listing for the dictionary " + str(dictionary_name) + " (" + bold_string(str(len(some_dictionary)) + " items)") + " *")
            print(title_string)
            print("")
 
            
            #print("print_dict function: type(some_dictionary) %s" %type(some_dictionary))
            if isinstance(some_dictionary, dict):
                for key in some_dictionary:
                    quoted_key, quoted_value = None, None

                    value = some_dictionary[ key ]
                    dictionary_name=str(dictionary_name)
                    quoted_key = "'" + str(key) + "'"
    
                    if isinstance(value, dict):
                        another_dictionary_detected = None
                        another_icon = str(dictionary_name) + "['" + str(key) + "'] "
                        len_of_value = len(value)
            
                        if len_of_value > 0:
                            print("")
                            print("")
                            output_statement2 = "         " + "{0: >20}".format("  %s dictionary has %s entry(s). The contents of this dictionary is:" %(str(another_icon), str(len_of_value)))
                            print((bold_string(output_statement2)))
                            print("")
                            if value != None:
                                print_inner_dictionary(value)
                            else:
                                print_info('print_dict','for key %s, the corresponding value is unset' %key)

                        else:
                            another_dictionary_detected = "  " + "{0: >31}".format(str(quoted_key)) + " => {}    (empty dictionary) "

                        end_statement = "                                                              ( end of output ) "
                        print(end_statement)
                        print("")
                        print("")

                    elif isinstance(value, list):
                        print("")
                        list_name =  str(key) 
                        another_icon =  " the value is of type list ['...'] "
                        len_of_value = len(value)
                        if len_of_value > 0:
                            list_detected = "       " + "{0: >20}".format(quoted_key) + "   =>    %s" %another_icon
                            #list_detected = "" + "{0: >20}".format(quoted_key) + " => %s" %another_icon
                            print((bold_string(list_detected)))
                            output_items = "" + "{0: >20}".format("") + "                The %s entries of this list is:" %str(len_of_value)
                            print(output_items)
                            #value_icon = "             " + "{0: >25}".format("  key    =>  ") + str(quoted_key) + "="
                            for element in value:
                                value_icon = "             " + "{0: >25}".format("                 =>    " + str(element))
                                print(value_icon)
                            #print "The attached list is:"
                            #print "value %s" %value
                            #print_list(value, 'value', 'padding')
                        else:
                            list_detected =  "  " + "{0: >25}".format(str(quoted_key)) + " => []    (empty list) "
                            print(list_detected)

                        end_statement = "                                    ( end of output for list %s ) " %list_name
                        print(end_statement)
                        print("")

                    elif isinstance(value, set):
                        set_name = str(key)
                        another_icon = " set ('" + str(key) + "') "
                        len_of_value = len(value)
                        if len_of_value > 0:
                            set_detected = "         " + "{0: >20}".format(quoted_key) + "    =>  %s" %another_icon
                            print((bold_string(set_detected)))
                            print("")
                            output_items = "" + "{0: >20}".format("") + "         The %s entries are:" %str(len_of_value)
                            print(output_items)
                            for element in value:
                                value_icon = "             " + "{0: >25}".format("  => " + str(element))
                                print(value_icon)
                        else:
                            set_detected =  "  " + "{0: >25}".format(str(quoted_key)) + " => []    (empty set) "
                            print(set_detected)
                        end_statement = "                                    ( end of output for set %s ) " %set_name
                        print(end_statement)
                        print("")

                    elif isinstance(value, str):
                        quoted_value = " '" + str(value) + "' "
                        output_statement = "   " + "{0: >25}".format(quoted_key) + "  =>  " + str(quoted_value)
                        print(output_statement)

                    else:
                        #print("else, print_dict-- begin")
                        #print("value %s" %value)
                        #print("type(value) %s" %type(value))
                        # Output of regulat key/value pairs.
                        quoted_value = None

                        if value != None:
                            if isinstance(value, int):
                                quoted_value = int(value)
                            elif isinstance(value, float):
                                quoted_value = float(value)

                            output_statement2 = "   " + "{0: >25}".format(quoted_key) + "  =>   " + str(quoted_value)
                            print(output_statement2)
                        #print("else, print_dict-- end")

                print("")
                end_string = "              *  end of " + str(dictionary_name) + " dictionary content listing"
                print((bold_string(end_string)))
                print("")

            else:
                print("")
                print_info('print_dict','ERROR: %s is not of type dictionary; it is type(some_dictionary) %s'  %(some_dictionary,type(some_dictionary)))
                print_info('print_dict','was the order reversed when invoking this function?')
                print_info('print_dict','format: print_dict(dictionary_variable, "dictionary_varible_name_in_quote")')
                print("")

    elif dictionary_name and not some_dictionary:
        print("")
        print_info('print_dict','unable to print the contents of dictionary %s since it is empty/None' %dictionary_name)
        print("")
    else:
        print_one_of_the_input_variables_not_set_error_and_exit('print_dict',some_dictionary,dictionary_name)


def print_inner_dictionary(some_dictionary):

    if some_dictionary:
        for key,value in list(some_dictionary.items()):
            quoted_key, quoted_value = None, None

            answer, key = is_it_an_integer_type(key)
            if answer == 'No':
                quoted_key = "'" + str(key) + "'"
            else:
                quoted_key = str(key)

            if isinstance(value, dict):
                another_dictionary_detected = None
                dictionary_name = quoted_key 
                len_of_value = len(value)
                inner_dictionary_icon = "dictionary with %s items: " %len_of_value
                if len_of_value > 0:
                    another_dictionary_detected = "                               " +  "{0: >38}".format(str(quoted_key)) + "  =>  %s" %inner_dictionary_icon
                    print((bold_string(another_dictionary_detected)))
                    print("")
                    for more_key, more_value in list(value.items()):
                        quoted_more_key, quoted_more_value = None, None

                        answer, more_key = is_it_an_integer_type(more_key)
                        if answer == 'No':
                            quoted_more_key = "'" + str(more_key) + "'"   
                        else:
                            quoted_more_key = str(more_key)

                        if more_value:
                            if not isinstance(more_value, dict) and not isinstance(more_value, list) and not isinstance(more_value, set):
                                answer, more_value = is_it_an_integer_type(more_value)
                                if answer == 'No':
                                    quoted_more_value = "'" + str(more_value)  + "'"
                                else:
                                    quoted_more_value = str(more_value) 
                            else:
                                quoted_more_value = more_value

                            content_of_inner_dict = "                                                           " +  "{0: >38}".format(quoted_more_key) + "  =>  " + str(more_value)
                            print(content_of_inner_dict)
                    print("") 
                else:
                    another_dictionary_detected = "                                            " +  "{0: >38}".format(quoted_key) + "  => {}    (empty dictionary) "
                    print(another_dictionary_detected)
            elif isinstance(value, list) or isinstance(value, set):
                inner_list_icon = None
                len_of_value = len(value)

                if isinstance(value, list):
                    inner_list_icon = " list has %s items: " %len_of_value
                else:
                    inner_list_icon = " set has %s items: " %len_of_value

                if len_of_value > 0:
                    another_list_detected = "                               " +  "{0: >38}".format(quoted_key) + "  =>  %s" %inner_list_icon
                    print((bold_string(another_list_detected)))
                    print("")
                    #output_items = "                            " + "{0: >38}".format("") + "         The %s entries are:" %str(len_of_value)
                    #print output_items
                    i = 1
                    for item in value:

                        if isinstance(item, dict):
                            len_of_item = len(item)
                            description = "                                         " + "{0: >38}".format("(" + str(i) + ")") + " dictionary with %s items: " %len_of_item
                            print(description)
                            if len_of_item > 0:
                                for key_item, value_item in list(item.items()): 
                                    content_of_inner_dict = "                                                                            " +  "{0: <38}".format(str(key_item) + "  =>  " + str(value_item))
                                    print(content_of_inner_dict)
                                print("")
                            else:
                                another_dictionary_detected = "                                                                            " +  "{0: <38}".format(str(item)) + " =>  {}    (empty dictionary) "
                                print(another_dictionary_detected)
                        elif isinstance(item, list) or isinstance(value, set):
                            description = None

                            len_of_item = len(item)
                            if isinstance(value, list):
                                description = "                                         " + "{0: >38}".format("(" + str(i) + ")") + " list with %s items: " %len_of_item
                            else:
                                description = "                                         " + "{0: >38}".format("(" + str(i) + ")") + " set with %s items: " %len_of_item

                            print(description)
                            #inner_list_icon = " list with %s items: " %len_of_item
                            if len_of_item > 0:
                                #another_list_detected = "                                       " +  "{0: >38}".format(str(item)) + " =>  %s" %inner_list_icon
                                #print another_list_detected
                                for an_item in item:
                                    content_of_inner_list = "                                       " +  "{0: >38}".format("         ") + " =>  " +  str(an_item)
                                    print(content_of_inner_list)
                            else:
                                another_list_detected = "                                       " +  "{0: >38}".format(str(item)) + " => []    (empty list) "
                                print(another_list_detected)
                        else:
                            output = "                                       " +  "{0: >38}".format(str(item))
                            print(output)
                        i = i + 1

                    #print output_items
                    #print_list(value, 'value', 'padding')

                else:
                    output = "                                  " +  "{0: >35}".format(quoted_key) + "  => []    (empty list) " 
                    print(output)
            else:
                if value:
                    answer, value = is_it_an_integer_type(value)
                    if answer == 'No':
                        quoted_value = "'" + str(value) + "'"
                    else:
                        quoted_value = str(value)
                    output_statement2 = "                                  " + "{0: >35}".format(str(quoted_key)) + "  => %s  " %quoted_value
                    print(output_statement2)
        print("")

    else:
        print_input_variable_not_set_error_and_exit('print_inner_dictionary','some_dictionary')

def print_dict_ordered(some_dictionary,some_dictionary_name):

    if some_dictionary and some_dictionary_name:
        print("")
        print_info('print_dict_ordered','dictionary name: %s' %bold_string(some_dictionary_name))
        print_info('print_dict_ordered','number of entries: %s' %bold_string(len(some_dictionary)))
        print("")
        ordered_key_list = sorted(some_dictionary)
        for ordered_key in ordered_key_list:
            value = some_dictionary[ordered_key]
            print(("  %s => %s " %(ordered_key,value)))

        print("")
        print("")

    else:
        print_one_of_the_input_variables_not_set_error_and_exit('run_bash_command',some_dictionary,some_dictionary_name)


def print_list(some_list, some_list_name,padding=None):
    number_of_elements, width = None, None

    if some_list:
        if isinstance(some_list,list):
            number_of_elements = len(some_list)
            print("")
            list_title = bold_string("              List %s has %s item(s), which are:" %(str(some_list_name), number_of_elements))
            print(list_title)


            if len(some_list) > 0:
                if int(number_of_elements) <= 6:
                    two_wide_row_print_stmt(some_list,padding=None)
                elif int(number_of_elements) > 6:
                    three_wide_row_print_stmt(some_list,padding=None)

        else:
            print("")
            print_info('print_list','input is not of type list-- perhaps the order was reversed')
            print_info('print_list','when invoking this function, it should be print_list(list_variable,"list_variable_name")')
            print("")
    else:
        print_one_of_the_input_variables_not_set_error_and_exit('print_list',some_list,some_list_name)


def two_wide_row_print_stmt(some_list,padding=None):
    # internal function called by the print_list function

    if some_list:
        list_size = len(some_list)
        highest_index_number = list_size - 1
        #print ""
        #print "  This list has %s members" %list_size
        print("")
        if highest_index_number == 0:
            if padding:
                print(("\t\t\t\t%s" %some_list[ 0 ]))
            else:
                print(("\t\t%s" %some_list[ 0 ]))
        else:
            for element in range(0,highest_index_number, 2):
                if element + 1 <=  highest_index_number:
                    if padding:        
                        print(("\t\t\t\t%s\t\t%s" %(some_list[ element ],  some_list[ element + 1 ] )))
                    else:
                        print(("\t\t%s\t\t%s" %(some_list[ element ],  some_list[ element + 1 ] )))
                else:
                    if padding:
                        print(("\t\t\t\t%s" %some_list[ element ]))
                    else:
                        print(("\t\t%s" %some_list[ element ]))

                # If highest_index_number is an odd number, print it also
                # since the for lop would end at the even number below it.
                if element + 2 == highest_index_number:
                    if padding:
                        print(("\t\t\t\t%s" %some_list[ element + 2 ]))
                    else:
                        print(("\t\t%s" %some_list[ element + 2 ]))
        print("")

    else:
        print_input_variable_not_set_error_and_exit('two_wide_row_print_stmt','some_list')


def three_wide_row_print_stmt(some_list,padding=None):
    # internal function called by the print_list function

    if some_list:
        list_size = len(some_list)
        highest_index_number = list_size - 1

        #print "   This list has %s members" %list_size
        #print ""

        for item in range(0,list_size, 3):
            if item + 3 <= list_size:
                if padding:
                    print(("\t\t\t\t\t%s\t\t%s\t\t%s" %(some_list[ item ], some_list[ item + 1 ], some_list[ item + 2 ] )))
                else:
                    print(("\t\t%s\t\t%s\t\t%s" %(some_list[ item ], some_list[ item + 1 ], some_list[ item + 2 ] )))
            elif item + 2 <= list_size:
                if padding:
                    print(("\t\t\t\t\t%s\t\t%s" %(some_list[ item ], some_list[ item + 1 ])))
                else:
                    print(("\t\t%s\t\t%s" %(some_list[ item ], some_list[ item + 1 ])))
            elif item + 1  <= list_size:
                if padding:
                    print(("\t\t\t\t%s" %some_list[ item ]))
                else:
                    print(("\t\t%s" %some_list[ item ]))
        print("")

    else:
        print_input_variable_not_set_error_and_exit('three_wide_row_print_stmt','some_list')

def delete_all_files_in_directory(some_dir):

    if some_dir:
        if os.path.exists(some_dir):
            # Avoid using shutil.rmtree:
            print_info('delete_all_file_in_directory','purging files in %s' %some_dir)

            for file_name in os.listdir(some_dir):
                if not file_name.startswith(".nfs"):
                    print_info('delete_all_file_in_directory','deleting file %s' %file_name)
                    os.remove(os.path.join(some_dir, file_name))
        else:
            print_path_does_not_exist_error_and_exit('delete_all_file_in_directory',some_dir)
    else:
        print_input_variable_not_set_error_and_exit('delete_all_files_in_directory','some_dir')

def parse_jamo_response_and_return_list(cmd):
    jamo_answer_list = []

    if cmd:
        #print_info('parse_jamo_response_and_return_list','running the command: %s' %cmd)
        output = run_bash_command(cmd)
        if output:
            jamo_answer_list = list(set(output.replace('None','').split('\n')))
            #print_info('parse_jamo_response_and_return_list','jamo_answer_list: %s' %jamo_answer_list)

            if jamo_answer_list:
                for jamo_answer_element in jamo_answer_list:
                    if jamo_answer_element == '':
                        jamo_answer_list.remove(jamo_answer_element)
            else:
                print_not_set_error_and_exit('parse_jamo_response_and_return_list','jamo_answer_list')
        else:
            print("")
            print_info('parse_jamo_response_and_return_list','no answer from jamo w/%s' %cmd)
            print("")
            
    else:
        print_input_variable_not_set_error_and_exit('parse_jamo_response_and_return_list','cmd_list')

    return jamo_answer_list

def add_readme_to_nmdc_source_tar_file_list(some_list, tarball_file_name, base_dir,file_format):
    ALT_BASE_DIR, centrifuge_readme_file, kraken2_readme_file, gottcha2_readme_file = None, None, None, None


    if some_list and tarball_file_name and base_dir and file_format:
        if os.environ.get('ALT_BASE_DIR'):
            ALT_BASE_DIR = os.environ.get('ALT_BASE_DIR')
        else:
            print_not_set_error_and_exit('add_readme_to_nmdc_source_tar_file_list','ALT_BASE_DIR')

        base_source_dir = str(ALT_BASE_DIR)

        if file_format == 'tsv_file_format':
            centrifuge_readme_file = str(ALT_BASE_DIR) + "/NMDC_readme/README.centrifuge.txt"
            kraken2_readme_file = str(ALT_BASE_DIR) + "/NMDC_readme/README.kraken2.txt"
            gottcha2_readme_file = str(ALT_BASE_DIR) + "/NMDC_readme/README.gottcha2.txt"
        elif file_format == 'csv_file_format':
            centrifuge_readme_file = str(ALT_BASE_DIR) + "/NMDC_readme_older_format/README.centrifuge.txt"
            kraken2_readme_file = str(ALT_BASE_DIR) + "/NMDC_readme_older_format/README.kraken2.txt"
            gottcha2_readme_file = str(ALT_BASE_DIR) + "/NMDC_readme_older_format/README.gottcha2.txt"
        else:
            print_info('add_readme_to_nmdc_source_tar_file_list','file_format is %s' %file_format)
            print_error('add_readme_to_nmdc_source_tar_file_list', 'much be either "tsv_file_format" or "csv_file_format"')

        if len(some_list) > 0:
            if 'centrifuge' in tarball_file_name:
                name = os.path.basename(centrifuge_readme_file)
                source_scratch_file_name = str(base_dir) + str(name)
                shutil.copy(centrifuge_readme_file,source_scratch_file_name)
                some_list.append(source_scratch_file_name)

            elif 'kraken2' in tarball_file_name:
                name = os.path.basename(kraken2_readme_file)
                source_scratch_file_name = str(base_dir) + str(name)
                shutil.copy(kraken2_readme_file,source_scratch_file_name)
                some_list.append(source_scratch_file_name)

            elif 'gottcha2' in tarball_file_name:
                name = os.path.basename(gottcha2_readme_file)
                source_scratch_file_name = str(base_dir) + str(name)
                shutil.copy(gottcha2_readme_file,source_scratch_file_name)
                some_list.append(source_scratch_file_name)

        else:
            print("")
            print_info('add_readme_to_nmdc_source_tar_file_list', 'the list is empty')
            print_info('add_readme_to_nmdc_source_tar_file_list', 'not adding the readme file to this list')
            print("")

    else:
        print_one_of_the_input_variables_not_set_error_and_exit('add_readme_to_nmdc_source_tar_file_list',some_list,tarball_file_name)
    
    return some_list

def untar_file(tarball_file_name):
    
    if tarball_file_name:
        directory_path = os.path.dirname(tarball_file_name)
        file_name = os.path.basename(tarball_file_name)
        if directory_path and file_name:
            if os.path.exists(tarball_file_name) and os.path.exists(directory_path):
                directory_list = []
                dir_name, base_dir, untarred_dirdctory_path = None, None, None

                # (1) save the cwd of the script that invoked this function.
                # (2) cd dest_directory
                # (3) untar if it has not been previously done so.
                # (4) cd back to directory where the script was located
                directory_location_of_function_invocation = os.getcwd()
                os.chdir(directory_path)
                tar = tarfile.open(file_name)
                tar.getmembers()
                for member in tar.getmembers():
                    # get each member file_name and check that all the files have a common subdirecotry.
                    # if so, we can then determine if the untarred directory already exists.

                    file_name_list = member.name.split('/')
                    if len(file_name_list) > 1:
                        dir_name = file_name_list[0] 
                        directory_list.append(dir_name)

                if len(list(set(directory_list))) == 1:
                    base_dir = list(set(directory_list))[0]

                if base_dir:
                    untarred_directory_path = str(directory_path) + "/" + str(base_dir)
                    if not os.path.exists(untarred_directory_path):
                        print_info('untar_file','extracting %s' %tarball_file_name)
                        tar.extractall()
                    else:
                        print_info('untar_file','the tarball was previously untarred to %s; no point untarring file' %untarred_directory_path)
                else:
                    tar.extractall()

                tar.close()

                print_info('untar_file','cleanup: deleting file %s' %tarball_file_name)
                os.remove(tarball_file_name)

                os.chdir(directory_location_of_function_invocation)

            else:
                print_path_does_not_exist_error_and_exit('untar_file',directory_path)
                print_path_does_not_exist_error_and_exit('untar_file',tarball_file_name)

        else:
            print_one_of_the_variables_not_set_error_and_exit('untar_file',directory,file_name)
    else:
        print_input_variable_not_set_error_and_exit('untar_file','tarball_file_name')

    
def copy_tarball_then_untar(full_path_src_tarball_file_name, destination_dir):

    if full_path_src_tarball_file_name and destination_dir:
        if os.path.exists(full_path_src_tarball_file_name):
            file_name = os.path.basename(full_path_src_tarball_file_name)
            destination_full_path_file_name = str(destination_dir) + "/" + str(file_name)

            if not os.path.exists(destination_full_path_file_name):
                shutil.copy(full_path_src_tarball_file_name, destination_full_path_file_name)
            else:
                print("")
                print_info('copy_tarball_then_untar','the file %s already exists, will not copy over it'  %destination_full_path_file_name)
                print("")

            if os.path.exists(destination_full_path_file_name):
                untar_file(destination_full_path_file_name)
            else:
                print_path_does_not_exist_error_and_exit('copy_tarball_then_untar',destination_full_path_file_name)
        else:
            print_path_does_not_exist_error_and_exit('copy_tarball_then_untar',full_path_src_tarball_file_name)
    else:
        print_one_of_the_input_variables_not_set_error_and_exit('copy_tarball_then_untar',full_path_src_tarball_file_name,destination_dir)


def create_tarball(some_file_name_dict, some_base_dir):
    # some_file_name_dict, a one entry dictionary containing:
    # some_file_name_dict[ tarball_name ] = [ list of source files to be tarred and gzipped ]
    #
    # some_base_dir = some scratch directory where the tarball is generated.

    generated_tarball_file_dict = {}

    if some_file_name_dict and some_base_dir:
        # the base_dir needs to have a trailing slash since the that will be used in a regular
        # expression to strip this full path when the tarball is created. add the traling '/'
        # and remove duplicate slashes if present.
        base_dir_path = some_base_dir + "/"
        base_dir = base_dir_path.replace('//', '/')

        # one of the / at the begining needs to be removed in order for this regex to work
        # example:
        # \/global\/cfs\/cdirs\/m342\/img\/nmdc_scratch_dir\/503125_159848\/ ->
        # \global\/projectb\/sandbox\/IMG\/img\/dataLoad6\/phajek\/tmp\/nmdc\/
        base_dir_delimited = base_dir.replace('/','\/').replace('\/global','\global')

        if len(some_file_name_dict) == 1:
            tarball_file_name = list(some_file_name_dict.keys())[0]
            source_scratch_file_list = some_file_name_dict[ tarball_file_name ]

            full_path_tarball_file_name = str(base_dir) + str(tarball_file_name)

            # pointless having the statements below since the dir containing these files would have been purged
            # but I am adding them as a reminder that this ontent need to be purged when this function is called.
            if os.path.exists(full_path_tarball_file_name):
                os.remove(full_path_tarball_file_name)

            if len(source_scratch_file_list) > 0:
                tar_filename = ' '.join(source_scratch_file_list)
                transform_string = "s/" + str(base_dir_delimited) + "//g"
                cmd = "cd %s && tar -cf %s %s --transform='%s'" %(base_dir, full_path_tarball_file_name, tar_filename, transform_string)
                print_info('create_tarball','generating tarball %s' %full_path_tarball_file_name)
                #print "  create_tarball: INFO: %s" %cmd
                run_bash_command(cmd)

                for file_name in source_scratch_file_list:
                    if not file_name.endswith(".tar") and file_name.endswith(".tar.gz"):
                        print_info('create_tarball','deleting file %s in scratch dir'  %file_name)
                        os.remove(file_name)
            else:
                print_not_set_error_and_exit('create_tarball','source_scratch_file_list')

            if os.path.exists(full_path_tarball_file_name):
                generated_tarball_file_dict[ tarball_file_name ] = full_path_tarball_file_name

        else:
            print_error('create_tarball','the input dictionary should one have length == 1')
    else:
        print_one_of_the_input_variables_not_set_error_and_exit('create_tarball',some_file_name_dict,some_base_dir)

    return generated_tarball_file_dict

def append_to_tarball(some_tarball_file_name, some_additional_file):

    if some_tarball_file_name and some_additional_file:
        if os.path.exists(some_tarball_file_name) and os.path.exists(some_additional_file):
            cmd = "tar -rf %s %s" %(some_tarball_file_name, some_additional_file)
            run_bash_command(cmd)
        else:
            print_path_does_not_exist_error_and_exit('append_to_tarball',some_tarball_file_name)
            print_path_does_not_exist_error_and_exit('append_to_tarball',some_additional_file)
    else:
        print_one_of_the_input_variables_not_set_error_and_exit('append_to_tarball',some_tarball_file_name,some_additional_file)

def gzip_file_name(file_name):

    if file_name:
        if os.path.exists(file_name):
            cmd = "gzip %s" %str(file_name)
            print("")
            print_info('gzip_file_name','running gzip %s' %file_name)
            print("")

            run_bash_command(cmd)

        else:
            print_path_does_not_exist_error_and_exit('gzip_file_name',file_name)
    else:
        print_input_variable_not_set_error_and_exit('gzip_file_name','file_name')

def get_file_name_list_from_tarball(file_name):
    file_name_list = []

    if file_name:
        if os.path.exists(file_name):
            if os.stat(file_name).st_size == 0:
                print("")
                print_info('check_the_files_in_the_tarballs','the file %s is a sparse file' %file_name)
                print("")
            else:
                tar = tarfile.open(file_name)
                for tar_file_member in tar.getmembers():
                    file_name_list.append(tar_file_member.name)             
        else:
            print_path_does_not_exist_error_and_exit('get_file_name_list_from_tarbal',file_name)
    else:
        print_input_variable_not_set_error_and_exit('get_file_name_list_from_tarball','file_name')

    return file_name_list

def is_it_an_mbin_oid_or_a_nmdc_key(some_string):
    taxon_oid_or_sp_gold_id, taxon_oid, bin_oid, nmdc_key, sp_gold_id, database = None, None, None, None, None, None
    mbin_list = []

    if some_string:
        from img_local_common_queries import check_and_validate_type_string
        taxon_oid_or_sp_gold_id, mbin_list, database = check_and_validate_type_string(some_string)
        if taxon_oid_or_sp_gold_id:
            if database:
                if database == 'mbin':
                    bin_oid = some_string
                    taxon_oid = taxon_oid_or_sp_gold_id
 
                    return bin_oid, taxon_oid, database

                elif database == 'nmdc':
                    nmdc_key = some_string
                    sp_gold_id = taxon_oid_or_sp_gold_id

                    return nmdc_key, sp_gold_id, database
    else:
        print_input_variable_not_set_error_and_exit('is_it_an_mbin_oid_or_a_nmdc_key','some_string')

    return None, None, None

def load_dictionary_from_file(file_name, separator, ordered=None):
    # loads file with two distinct values and loads it into the dictionary
    # some _dict. The file value in the line will represent the key of the dictionary and
    # the second value in the line will be the value of this dictionary.
    #
    # I added ordered=None input parameter so we can recieve an ordered dictionary if the
    # filed is set. -phajek
    some_dict = None

    if file_name and separator:
        if ordered == None:
            some_dict = {}
        else:
            some_dict = OrderedDict()

        if os.path.exists(file_name):
            with open(file_name, "r") as in_fd:
                for line in in_fd.readlines():
                    if line.startswith('#'):
                        # skip header line
                        next(in_fd)
                    else:
                        key_pair = line.rstrip().split(separator)

                        if len(key_pair) == 2:
                            some_dict[ convert_some_value_into_more_usable_format(key_pair[0]) ] = convert_some_value_into_more_usable_format( key_pair[1] )

                        elif len(key_pair) > 2:
                            first_element = key_pair[0]
                            del key_pair[0]
                            the_remaining_string = separator.join(key_pair)
                            some_dict[ first_element ] = the_remaining_string
        else:
            print_path_does_not_exist_error_and_exit('load_dictionary_from_file',file_name)
    else:
        print_one_of_the_input_variables_not_set_error_and_exit('load_dictionary_from_file',file_name,separator)

    return some_dict

def load_taxon_stats_from_merfs_filesystem(accounting):
    taxon_stats_dict = {}

    if accounting:
        if accounting.get('dataLoad_path'):
            dataLoad_path = accounting[ 'dataLoad_path' ]
        else:
            print_not_set_error_and_exit('load_taxon_stats_from_merfs_filesystem','dataLoad_path')

        if accounting.get('submission_id'):
            submission_id = accounting[ 'submission_id' ]

        if  accounting.get('taxon_oid'):
             taxon_oid = accounting[ 'taxon_oid' ]
        else:
            print_not_set_error_and_exit('load_taxon_stats_from_merfs_filesystem','taxon_oid')

        base_dir = str(dataLoad_path) + '/web-data/mer.fs/' + str(taxon_oid) + '/assembled'
        if os.path.exists(base_dir):
            accounting[ 'base_dir' ] = base_dir
            merfs_taxon_stats_file = str(base_dir) + '/taxon_stats.txt'
            if os.path.exists(merfs_taxon_stats_file):
                taxon_stats_dict = load_dictionary_from_file(merfs_taxon_stats_file, '\t')
            else:
                print_path_does_not_exist_error_and_exit('load_taxon_stats_from_merfs_filesystem',merfs_taxon_stats_file)

            if taxon_stats_dict:
                if len(taxon_stats_dict) > 0:
                    accounting[ 'taxon_stats' ] = taxon_stats_dict

                elif len(taxon_stats_dict) == 0:
                    print_error('load_taxon_stats_from_merfs_filesystem','the taxon_stats dictionary is empty')
            else:
                print_not_set_error_and_exit('load_taxon_stats_from_merfs_filesystem','taxon_stats')
        else:
            print_path_does_not_exist_error_and_exit('load_taxon_stats_from_merfs_filesystem',base_dir)
    else:
        print_input_variable_not_set_error_and_exit('load_taxon_stats_from_merfs_filesystem','accounting')

    return accounting

def get_file_line_count(some_file):
    c_generator = None

    if some_file:
        if os.path.exists(some_file):
            f = open( some_file, 'rb')
            c_generator = _count_generator(f.raw.read)

            return sum( buf.count(b'\n') for buf in  c_generator)

        else:
            print_path_does_not_exist_error_and_exit('get_file_line_count',some_file)
    else:
        print_input_variable_not_set_error_and_exit('get_file_line_count','some_file')

    return line_count

def _count_generator(reader):
    # used with get_file_line_count function
    
    if reader:
        b = reader(1024 * 1024)

        while b:
            yield b
            b = reader(1024 * 1024)
    else:
        print_input_variable_not_set_error_and_exit('_count_generator','reader')


def read_first_line_of_file(some_file):
    # returns a list of the the first line

    first_line = None
    first_line_list = []

    if some_file:
        if os.path.exists(some_file):
            with open(some_file, "r") as in_fd:
                for line in in_fd.readlines():
                    line_content = line.replace("\n", "")
                    line_list = line_content.split('\t')
                    for some_line in line_list:
                        first_line_list.append(convert_some_value_into_more_usable_format(some_line))
                        return first_line_list
        else:
            print_path_does_not_exist_error_and_exit('read_first_line_of_file', some_file)
    else:
        print_input_variable_not_set_error_and_exit('read_first_line_of_file','some_file')


def read_from_file_to_get_a_list(some_file):
    # Often we have a file containing a list of submission_ids or taxon_oids 
    # with each line containing a submission_id, taxon_oid, etc. The functions reads the file and returns 
    # a list of non-null content (new line is the list delimiter.
    some_list = []

    if some_file:
        if os.path.exists(some_file):
            # os.stat empty file = 0, populated_file = 1
            result = os.stat(some_file)
            if result != 0:
                with open(some_file, "r") as in_fd:
                    for line in in_fd.readlines():
                        line_content = line.replace("\n", "")
                        line_correctly_typed = convert_some_value_into_more_usable_format(line_content)
                        some_list.append(line_correctly_typed)
            else:
                print_info('read_from_file_to_get_a_list','the file %s is empty' %some_file)
                return []
        else:
            print_path_does_not_exist_error_and_exit('read_from_file_to_get_a_list',some_file)
    else:
        print_input_variable_not_set_error_and_exit('read_from_file_to_get_a_list','some_dictionary')

    return  some_list


def copy_file_to_dna_volume(some_infile,some_outfile_dir):
    # The dna volume (which the UI sources) is only writing by
    # root hosts like dtn04.nersc.gov and not w/most hosts including the NERSCs VMs.
    file_name, dna_file_name, copy_command = None, None, None

    if some_infile and some_outfile_dir:
        dtn04_host = "dtn04"

        if not os.path.exists(some_infile):
            print_path_does_not_exist_error_and_exit('copy_file_to_dna_volume',some_infile)

        else:
            file_name = os.path.basename(some_infile)
            dna_file_name = str(some_outfile_dir) + '/' + str(file_name)
            copy_command = 'ssh -q ' + str(dtn04_host) + ' cp ' + str(some_infile) + ' ' + str(some_outfile_dir)

            if os.path.exists(dna_file_name):
                print("")
                print_info('copy_file_to_dna_volume','deleting the previously generated file %s' %dna_file_name)
                print("")

                cmd = "ssh -q dtn04 rm %s" %dna_file_name

                print("")
                print_info('copy_file_to_dna_volume','running the command %s' %cmd)
                print("")

                run_bash_command(cmd)
                if os.path.exists(dna_file_name):
                    while os.path.exists(dna_file_name):
                        print("")
                        print_info('copy_file_to_dna_volume','waiting for the removal of %s to finish on dtn04 and be visible on all systems' %dna_file_name)
                        print("")
                        time.sleep(2)

            if os.path.exists(some_outfile_dir):
                if not os.path.exists(dna_file_name):
                    print("")
                    print_info('copy_file_to_dna_volume','now copying %s to %s' %(dna_file_name,some_outfile_dir))
                    print("")

                    run_bash_command(copy_command)
                    sleep_x_times_waiting_for_file_presence_on_the_file_system(dna_file_name,'wait_til_done')


            elif not os.path.exists(some_outfile_dir):
                cmd = "ssh -q dtn04 mkdir %s" %some_outfile_dir

                print_info('copy_file_to_dna_volume','running the command %s' %cmd)
                print("")

                run_bash_command(copy_command)
                sleep_x_times_waiting_for_file_presence_on_the_file_system(some_outfile_dir,'wait_til_done')


                if not os.path.exists(dna_file_name):

                    print_info('copy_file_to_dna_volume','now copying %s to %s' %(dna_file_name,some_outfile_dir))
                    print("")

                    run_bash_command(copy_command)
                    sleep_x_times_waiting_for_file_presence_on_the_file_system(dna_file_name,'wait_til_done')

    else:
        print_one_of_the_input_variables_not_set_error_and_exit('copy_file_to_dna_volume',some_infile,some_outfile_dir)
     

def get_todays_date(date_format = None):
    todays_date = None

    if date_format:
        #print "  get_todays_date: INFO: used date format %s" %date_format
        todays_date = datetime.date.today().strftime(date_format)
    else:    
        todays_date = datetime.datetime.now().strftime("%d-%m-%y")

    return todays_date

def prep_variable_for_dml_construction(some_var_name, some_var, dml_type):
    # example: prep_variable_for_dml_construction(its_spid, 'its_spid', 'update')
    #           where its_spid = 1228746
    #
    # output:
    #           its_spid_sql = 'its_spid='
    #           its_spid_sql2 = '1228746,'
    # 
    # if a variable is None, set its_spid_sql = '' and its_spid_sql2 = ''
    #
    # The purpose is to dynamically generate the SQL string statement so we
    # don't need to generate all the different insert/update SQL based if the variable is populated or not.
    # if it is not populated, its_spid_sql + its_spid_sql2 =''
    # and it would not appear during the construction of the DML statement.
    some_variable_sql, some_variable_sql2 = None, None
    if some_var_name != None and dml_type != None:
        if dml_type == 'insert':
            if some_var != None:
                some_variable_sql = str(some_var_name) + ","

                cleaned_var = convert_some_value_into_more_usable_format(some_var)

                if isinstance(cleaned_var, float) or isinstance(cleaned_var, int):
                    some_variable_sql2 = str(some_var) + ","

                if isinstance(cleaned_var, str):
                    some_variable_sql2 = "'" + str(some_var) + "',"

                return some_variable_sql, some_variable_sql2

            else:
                some_variable_sql = ''
                some_variable_sql2 = ''

                return some_variable_sql, some_variable_sql2


        elif dml_type == 'update':
            if some_var != None:
                some_variable_sql = str(some_var_name) + "= "

                cleaned_var = convert_some_value_into_more_usable_format(some_var)
                if isinstance(cleaned_var, float) or isinstance(cleaned_var, int):
                    some_variable_sql2 = str(some_var) + ", "
                elif isinstance(cleaned_var, str):
                    some_variable_sql2 = "'" + str(some_var) + "',"
                else:
                    print_error('prep_variable_for_dml_construction', 'variable type %s is not supported here' %type(cleaned_var))

                return some_variable_sql, some_variable_sql2

            else:
                some_variable_sql = ''
                some_variable_sql2 = ''

                return some_variable_sql, some_variable_sql2
    else:
        print_one_of_the_input_variables_not_set_error_and_exit('prep_variable_for_dml_construction',some_variable_name,  dml_type)


def prep_date_variable_for_dml_construction(some_variable_name, some_variable, date_format, dml_type):
    some_variable_sql, some_variable_sql2 = None, None

    if some_variable_name and date_format and dml_type:
        if dml_type == 'insert':
            if some_variable:
                some_variable_sql = str(some_variable_name) + ", "
                some_variable_sql2 = "to_date('" + str(some_variable) + "', '" + str(date_format) + "'), "
            else:
                some_variable_sql = ''
                some_variable_sql2 = ''

            return some_variable_sql, some_variable_sql2

        elif dml_type == 'update':
            if some_variable:
                some_variable_sql = str(some_variable_name) + "= "
                some_variable_sql2 =  "to_date('" + str(some_variable) + "', '" + str(date_format) + "'), "
            else:
                some_variable_sql = ''
                some_variable_sql2 = ''

            return some_variable_sql, some_variable_sql2

    else:
        print_one_of_the_input_variables_not_set_error_and_exit('prep_date_variable_for_dml_construction',some_variable_name, dml_type)


def fetch_gold_id_from_string(some_string):
    gold_id = None
    if some_string:
        if some_string.startswith('Ga'):
            substring = convert_some_value_into_more_usable_format(some_string.replace('Ga',''))
            if isinstance(substring, int):
                return some_string

        if '.' in some_string:
            string_list = some_string.split('.')
            for sub_string in string_list:
                match = re.match(r'^(Ga[0-9]+)', sub_string)
                if match:
                    gold_id = match.group(1)
                    if gold_id:
                        return gold_id

        elif '_' in some_string:
            string_list = some_string.split('_')
            for sub_string in string_list:
                match = re.match(r'^(Ga[0-9]+)', sub_string)
                if match:
                    gold_id = match.group(1)
                    if gold_id:
                        return gold_id

    else:
        print_input_variable_not_set_error_and_exit('fetch_gold_id_from_string','some_string')


def fetch_sp_gold_id_from_string(some_string):
    sp_gold_id = None

    if some_string:
        if '.' in some_string:
            string_list = some_string.split('.')
            for sub_string in string_list:
                match = re.match(r'^(Gp[0-9]+)', sub_string)
                if match:
                    sp_gold_id = match.group(1)
                    if sp_gold_id:
                        return sp_gold_id


        elif '_' in some_string:
            string_list = some_string.split('_')
            for sub_string in string_list:
                match = re.match(r'^(Ga[0-9]+)', sub_string)
                if match:
                    sp_gold_id = match.group(1)
                    if sp_gold_id:

                        return sp_gold_id

    else:
        print_input_variable_not_set_error_and_exit('fetch_sp_gold_id_from_string','some_string')


def get_key_list(dict):
    return list(dict.keys())

def get_values_list(dict):
    return list(dict.values())

def check_file_name_path(some_file_name):
    # This function is to deterimine if the input value is a file_name/path or not.
    # Is so, sets output variables dir_path, file_name.
    #
    # If it is not a valid file, dir_path = 'invaid_path'
    # while passing the input variable as the file_name.
    #
    # Type of valid file_name as an input.
    # 3 types of file information. Using this python script as an example:
    #
    #   (1)             img_local_utilities.py
    #   (2)             ./img_local_utilities.py
    #   (3)             str(ALT_BASE_DIR) + /dataLoad/custom_python_modules/utilities/img_local_utilities.py
    #
    # In each of the 3 cases, output will be in 2 variables-- directory_path and file_name
    # where the full path would be:
    # 
    # full_file_name_path = str( directory_path ) + '/' + str(file_name)
    # 
    current_working_directory = os.getcwd()

    if some_file_name:
        if os.path.exists(some_file_name):
            file_name = os.path.basename(some_file_name)
            dir_path = os.path.dirname(some_file_name)
    
            if file_name and dir_path:
                full_file_name_path = str(dir_path) + '/' + str(file_name)

                return dir_path, file_name

            else:
                return 'invaid_path', 'invalid_file'
        else:
            return 'invaid_path', 'invalid_file'
    else:
        print_input_variable_not_set_error_and_exit('check_file_name_path','some_file_name')


def parse_supplied_arguments(user_supplied_argument_list, min_num_of_args_required, error_message ='<file_name> <submission_id|taxon_oid>'):
    # (1) If you are expecting some string value to be provided (such as a file name) along with either a submission_id
    #     or a taxon_oid, you can use this function to:
    #           (a) validate that the integer argument provided is either a submission_id or taxon_oid
    #               (i) Provide the corresponding taxon_oid or submission_id for that integer argument (if it exists)
    #               (ii) Throw an error and exit if some sort of logical rule has been violated, e.g., the provided integer 
    #                    is a taxon_oid, but there is not a submission_id mapped to it. 
    #
    # (2) Order of arguments for this function is irrelevant. The function will always grab the integer argument first, 
    #     identify whether it is a taxon_oid or a submission_id, and then grab its corresponding value before interacting
    #     with the string.
    #
    # (3) To call this function, you need to have the following:
    #           (a) The sys.argv values in a list variable.
    #           (b) min_num_of_args_required set to either 2 or 3.
    #
    # (4) To invoke this function with only one argument, do this:
    #           min_num_of_args_required = 2
    #           submission_id, taxon_oid, database = parse_supplied_arguments(user_supplied_argument_list, min_num_of_args_required)
    #
    # (5) To invoke this function with two arguments, do this:
    #           min_num_of_args_required = 3
    #           submission_id, taxon_oid, database, file_name = parse_supplied_arguments(user_supplied_argument_list,min_num_of_args_required)
    #
    # (6) To parse three arguments, you can override the default message of '<file_name> <submission_id|taxon_oid>' to
    #           something else, such as:
    #               '<submission_id|taxon_oid> <file_name>' 
    #           by adding it to the last argument, like so:
    #               submission_id, taxon_oid, database, file_name = parse_supplied_arguments(user_supplied_argument_list,min_num_of_args_required,'<submission_id|taxon_oid> <file_name>')
    #

    usage_message, some_string, is_number, submission_id, taxon_oid, database = None, None, None, None, None, None

    if user_supplied_argument_list and min_num_of_args_required:
        number_of_arguments_provided = len(user_supplied_argument_list) + 1
        script_name = str(user_supplied_argument_list[0])

        # Take care of the times where the user did not provide enough arguments to meet their own criteria. 
        if number_of_arguments_provided <= min_num_of_args_required:
            if min_num_of_args_required == 2:
                usage_message = ' * Usage: ' + str(script_name) + ' ' + '<submission_id|taxon_oid>'
                print(usage_message)
                sys.exit(1)

            elif min_num_of_args_required == 3:
                if  error_message == '<file_name> <submission_id|taxon_oid>':
                    usage_message = ' * Usage: ' + str(script_name) + ' ' + '<file_name> <submission_id|taxon_oid>'
                    print(usage_message)
                    sys.exit(1)

                else:
                    # If another error message was explicitly set on how to provided the parameters and in what order
                    usage_message = ' * Usage: ' + str(script_name) + ' ' + str(error_message)
                    print(usage_message)
                    sys.exit(1)

        for i in range(len(user_supplied_argument_list)):
            if i > 0:
                some_element = user_supplied_argument_list[i]
                answer, some_element = is_it_an_integer_type(some_element)
                if answer == 'Yes':
                    # If taxon_oid is valid but there is no submission_id, taxon_oid_or_submission_id will set taxon_oid = None
                    submission_id, taxon_oid, database = taxon_oid_or_submission_id(some_element)

                    if submission_id:
                        if submission_id == int(some_element):
                            if taxon_oid:
                                print("")
                                #print "parse_supplied_arguments: with a corresponding taxon_oid of %s" %taxon_oid
                                #print ""

                            elif not taxon_oid:
                                if database != 'IMG_MER_RNASEQ' and database != 'IMG_ER_RNASEQ':
                                    print_error('parse_supplied_arguments','something is wrong, it should have a taxon_oid mapped to submission_id %s but it DOES NOT.' %submission_id)
                    else:
                        if not taxon_oid and not submission_id:
                            print_error('parse_supplied_arguments','%s is neither a valid submission_id or taxon_oid' %user_supplied_argument_list[i])
                else:
                    if some_string:
                        print_error('parse_supplied_arguments','more than one non-integer was supplied as command line arguments to %s' %script_name)
                    else:
                        some_string = user_supplied_argument_list[i]
        if some_string:
            if len(user_supplied_argument_list) == 2:
                some_string = user_supplied_argument_list[1]
                return some_string
        else:
            return submission_id, taxon_oid, database

    else:
        print_one_of_the_variables_not_set_error_and_exit('parse_supplied_arguments',user_supplied_argument_list,min_num_of_args_required)


def is_it_an_integer_type(some_possible_value):
    # returns <'Yes'|'No'> and if type 'int', the variable is set int(some_value) or if type 'str', the variable is 
    # set str(some_value) 
    # is_it_an_integer_type should be called 'is_it_an_number_type'

    # for some stupid reason, in an if block, i.e. 'if some_value:', python treats None and 0 ( 0.0 ) similarly and will
    # not enter the if statement. 
    answer = None

    if some_possible_value != None:
        some_value = some_possible_value
        #print("is_it_an_integer_type function")
        #print("type(some_value) %s" %type(some_value))
        if some_value != None:
            if not isinstance(some_value, dict) and not isinstance(some_value, list) and not isinstance(some_value, set):
                try:
                    if isinstance(int(some_value), int):
                        return 'Yes', int(some_possible_value)
                except:
                    try:
                        if isinstance(float(some_value), float):
                            return 'No', float(some_possible_value)

                    except:
                        #print("checking str scenario")
                        if isinstance(some_value, str):
                            if some_value.isdigit():
                                return 'Yes', int(some_possible_value)
                            else:
                                return 'No', some_possible_value
                        else:
                            return 'No', some_possible_value
            else:
                return 'No', some_possible_value
        else:
            print_not_set_error_and_exit('is_it_an_integer_type','some_value')

    else:
        print_input_variable_not_set_error_and_exit('is_it_an_integer_type','some_possible_value')


def remove_dups(list_of_something):
    unique_list = []

    if list_of_something:
        for element in list_of_something:
            if element not in unique_list:
                unique_list.append(element)
    else:
        print_input_variable_not_set_error_and_exit('remove_dups','list_of_something')

    return unique_list

def create_dictionary_from_two_lists(some_list1, some_list2):
    # the first list above (some_list1) will be the keys, 
    #  and the second list (some_list2) will be the values.
    some_dictionary_tmp, some_dictionary = {}, {}

    if some_list1 and some_list2:
        if len(some_list1) == len(some_list2):
            some_dictionary = dict(list(zip(some_list1, some_list2)))

        else:
            print_info('create_dictionary_from_two_lists','the two lists do not have the same number of entries')
            print_error('create_dictionary_from_two_lists','i.e. len(some_list1): %s, len(some_list2): %s' %(some_list1, some_list2))
    else:
        print_one_of_the_input_variables_not_set_error_and_exit('create_dictionary_from_two_lists',some_list1,some_list2)

    return some_dictionary

def load_json_config(accounting):
    dataLoad_path, taxon_oid = None, None

    if accounting:
        if accounting.get('dataLoad_path'):
            dataLoad_path = accounting.get('dataLoad_path')
        else:
            print_not_set_error_and_exit('load_json_config','dataLoad_path')

        if accounting.get('taxon_oid'):
            taxon_oid = accounting.get('taxon_oid')
        else:
            print_not_set_error_and_exit('load_json_config','taxon_oid')

        hostname = socket.gethostname()
        metagenomics_pipeline_status_file = str(dataLoad_path) + 'IMG_M_ER_pipeline_current_running_on_' + str(hostname)
        merfs_start_oid_file = str(dataLoad_path) + "/data/merfs/merfs_start_oid.txt"

        if os.path.exists(metagenomics_pipeline_status_file):
            # When used inline with the IMG/M ER Pipeline, the taxon_oid file created should match
            # the taxon_oid.  The safegaurds are needed due to one time when the taxon_oid was read
            # from the json file generated by the pipeline for another Submission ID.

            if os.path.exists(merfs_start_oid_file):
                stored_taxon_oid = None

                with open(taxon_oid_file, "r") as in_fd:
                    for line in in_fd.readlines():
                        stored_taxon_oid = int(line)

                if int(stored_taxon_oid) == int(taxon_oid):
                    print_info('load_json_config','loading the <submission_id>.json file into the accounting dict')

                    config_file_json = str(dataLoad_path) + "/data/merfs/submissions/" + str(sub_id) + ".json"
                    if os.path.exists(config_file_json):
                        with open (config_file_json) as f:
                            accounting = json.load(f)

                    # No sure if accounting[ 'dataLoad_path' ] is overwritten
                    accounting[ 'dataLoad_path' ] = dataLoad_path
                else:
                    print_error('load_json_config','unable to load the json file since the stored taxon_oid %s does not match the taxon_oid %s' %(stored_taxon_oid,taxon_oid))
            else:
                print_path_does_not_exist_error_and_exit('load_json_config','merfs_start_oid_file')
    else:
        print_input_variable_not_set_error_and_exit('load_json_config','accounting')

    return accounting


def load_rnaseq_expression_listfile_to_dict(file_name):
    # This script opens up the file given and attempts to create a dictionary with its_spid as 
    # keys and expression file(s) as values. This has been updated so that now each its_spid 
    # can hold more than one expression file, as we want to include intergenic expression data now. 
    
    rnaseq_expression_dict = {}

    if file_name:
        if os.path.exists(file_name):
            with open(file_name) as rnaseq_expression_lists:
                line_reader = csv.reader(rnaseq_expression_lists, delimiter='\t')
                next(line_reader)
                for line in line_reader:
                    its_spid, rnaseq_expression_file = None, None

                    its_spid = int(line[1])
                    rnaseq_expression_file = line[3]
                    if its_spid not in rnaseq_expression_dict:
                        rnaseq_expression_dict[ its_spid ] = [rnaseq_expression_file]
                    else:
                        rnaseq_expression_dict[ its_spid ].append(rnaseq_expression_file)
    else:
        print_input_variable_not_set_error_and_exit('load_rnaseq_expression_listfile_to_dict','file_name')

    return rnaseq_expression_dict


def taxon_oid_or_submission_id(some_number):
    # This script will check to see if a given number is a submission_id or a taxon_oid. It will then 
    # return the submission_id, taxon_oid, and database. If it is not found to be a submission_id or 
    # or taxon_oid in the submission_table, the script will check the img_core_v400.taxon table since 
    # some older submissions can be present there and not in the submission table. 
    
    submission_id, database, taxon_oid, sql = None, None, None, None

    if some_number:
        dsn_tns = cx_Oracle.makedsn('oracle-vm-1.jgi.lbl.gov', '1521', service_name='IMGIPRDPDB')
        connection = cx_Oracle.connect(user=r'imgsg_dev', password='Tuesday', dsn=dsn_tns)
        cursor = connection.cursor()
        if connection and cursor:
            answer, some_number = is_it_an_integer_type(some_number)

            # If the value given to the function is not an integer, then the script will not 
            # return anything and it will probably break whatever function is using it without 
            # giving a good explanation. I don't know if that's ever going to happen, but
            # it would probably be good to have some sort of contingency.
            if answer == 'Yes':
                sql = "select submission_id from submission where submission_id = %s" %int(some_number)
                cursor.execute(sql)
                if cursor.fetchone():
                    cursor.execute(sql)
                    # Also another really annoying issue:
                    # https://stackoverflow.com/questions/48143659/cursor-fetchone-returns-none-even-though-a-value-exists
                    # * cursor.fetchone() when called, does not retain the value so you have to perform another 
                    # * cursor.execute(sql) command.
                    #       yeah it blows.
                    #           Geez, Patrick. Bad day, huh? 
                    submission_id = cursor.fetchone()[0]
                    sql = "select database from submission where submission_id = %s" %int(submission_id)
                    cursor.execute(sql)
                    if cursor.fetchone():
                        cursor.execute(sql)
                        database = cursor.fetchone()[0]

                        if database != 'IMG_MER_RNASEQ' and  database != 'IMG_ER_RNASEQ':
                            sql = "select img_taxon_oid from submission where submission_id = %s" %int(submission_id)
                            cursor.execute(sql)

                            if  cursor.fetchone():
                                cursor.execute(sql)
                                taxon_oid = cursor.fetchone()[0]
                                # All three should have non-null values. 
                                return submission_id, taxon_oid, database

                            else:
                                print("")
                                print_info('taxon_oid_or_submission_id','since %s is a submission_id and from database %s, it should have a taxon_oid mapped to it' %(some_number, database))
                                print("")
                                # Returns submission_id, None, database
                                return submission_id, taxon_oid, database
                        else:
                            # Since it is an RNASeq submission_id, it will not have a taxon_oid mapped to it.
                            # Returns submission_id, None, database
                            return submission_id, taxon_oid, database
                    else:
                        print("")
                        print_info('taxon_oid_or_submission_id','for some_value %s, it has a submission_id %s but no database value mapped to that submission_id' %(some_number,submission_id))
                        print("")
                        # Returns submission_id, None, None
                        return submission_id, taxon_oid, database
                else:
                    # So far, this number is not a submission_id. Checking for it being a taxon_oid value.
                    sql  = "select img_taxon_oid from submission where img_taxon_oid = %s" %int(some_number)
                    cursor.execute(sql)
                    if cursor.fetchone():
                        cursor.execute(sql)
                        taxon_oid = cursor.fetchone()[0]
                        sql = "select submission_id from submission where img_taxon_oid =  %s" %taxon_oid
                        cursor.execute(sql)

                        if cursor.fetchone():
                            cursor.execute(sql)
                            submission_id = cursor.fetchone()[0]
                            sql = "select database from submission where submission_id = %s" %int(submission_id)
                            cursor.execute(sql)

                            if cursor.fetchone():
                                cursor.execute(sql)
                                database = cursor.fetchone()[0]
                                # All three should have non-null values. 
                                return submission_id, taxon_oid, database
                            else:
                                print("")
                                print_info('taxon_oid_or_submission_id','%s was determined to the taxon_oid, %s is the submission_id-- yet, there is no database parameter set in the submission table')
                                print("")
                                # Returns submission_id, taxon_oid, None
                                return submission_id, taxon_oid, database
                        else:
                            print("")
                            print_info('taxon_oid_or_submission_id','%s was determined to the taxon_oid, but there is no submission_id mapped to it' %some_number)
                            print("")
                            # Returns None, taxon_oid, None 
                            return submission_id, taxon_oid, database
                    else:
                        # For older taxon_oids, there may be no entry in the submission table, but still have 
                        # an entry in the img_core_v400.taxon table. That is what we are going to check. 
                        
                        sql = "select taxon_oid from img_core_v400.taxon where taxon_oid = %s" %some_number
                        cursor.execute(sql)
                        taxon_oid = cursor.fetchone()[0]
                        if taxon_oid:
                            sql = "select genome_type from img_core_v400.taxon where taxon_oid = %s" %taxon_oid
                            cursor.execute(sql)
                            genome_type = cursor.fetchone()[0]
                            if genome_type:
                                sql = "select submission_id from img_core_v400.taxon where taxon_oid = %s" %taxon_oid
                                cursor.execute(sql)
                                submission_id = cursor.fetchone()[0]
                                if genome_type == 'isolate':
                                    # Will return all non-null values.
                                    return submission_id, taxon_oid, 'IMG ER'
                                elif genome_type == 'metagenome':
                                    # Will return all non-null values.
                                    return submission_id, taxon_oid, 'IMG/M ER'
                                else:
                                    # Returns submission_id, taxon_oid, None.
                                    return submission_id, taxon_oid, None
                            else:
                                # Something weird has happened and we're not sure if it's possible. The message is here, but
                                # we do not have a plan to deal with the situation. 
                                print_info('taxon_oid_or_submission_id','odd, the taxon_oid exists in the taxon table yet the genome_type value is not set')
                        else:
                            # Check to see if it's a submission_id in the img_core_v400.taxon table. 
                            sql = "select taxon_oid from img_core_v400.taxon where submission_id = %s" %some_number
                            cursor.execute(sql)
                            taxon_oid = cursor.fetchone()[0]
                            if taxon_oid:
                                sql = "select genome_type from img_core_v400.taxon where taxon_oid = %s" %taxon_oid
                                cursor.execute(sql)
                                genome_type = cursor.fetchone()[0]
                                if genome_type:
                                    if genome_type == 'isolate':
                                        # Should return all non-null values.
                                        return submission_id, taxon_oid, 'IMG ER'
                                    elif genome_type == 'metagenome':
                                        # Should return all non-null values.
                                        return submission_id, taxon_oid, 'IMG/M ER'
                                else:
                                    # Returns submission_id, taxon_oid, None
                                    return submission_id, taxon_oid, None
                            else:
                                return None, None, None

        else:
            print_unable_to_connect_to_db('taxon_oid_or_submission_id','imgsg_dev')
    else:
        print_input_variable_not_set_error_and_exit('taxon_oid_or_submission_id','some_number')

def list_to_string(one_list):
    # should beusing the ohter function above.
    one_long_string = None

    if one_list:
        if isinstance(one_list, list):
            one_long_string = ','.join([str(element) for element in one_list])
        elif isinstance(one_list, str):
            one_long_string  = one_list
    else:
        print_input_variable_not_set_error_and_exit('list_to_string','one_list')

    return one_long_string

def clean_int_list(some_list):
    unique_item_list, cleaned_list = [], []
    int_or_str = ''
    
    if some_list:
        for element in some_list:
            if isinstance(element, str) or isinstance(element, int):
                if str(element).isdigit():
                    int_or_str = 'integer'
                    cleaned_list.append(int(element))

        if int_or_str == 'integer':
            unique_item_list = _remove_dups(cleaned_list)
        else:
            unique_item_list = _remove_dups(some_list)
    else:
        print_input_variable_not_set_error_and_exit('clean_int_list','some_list')

    return unique_item_list

def string_to_list(some_string):
    # for string containing ','
    if some_string:
        some_string  = str(some_string).rstrip(',')
        some_string  = str(some_string).lstrip(',')

        some_list = str(some_string).split(',')
        cleaned_list, unique_item_list = [], []
        int_or_str = ''

        for element in some_list:
            if element.isdigit():
                #print "int .. "
                int_or_str = 'integer'
                cleaned_list.append(int(element))

        if int_or_str == 'integer':
            unique_item_list = _remove_dups(cleaned_list)
        else:
            unique_item_list = _remove_dups(some_list)
    else:
        print_input_variable_not_set_error_and_exit('string_to_list','some_string')

    return unique_item_list


def sleep_x_times_waiting_for_file_presence_on_the_file_system(file_name,sleep_count):

    if file_name and sleep_count:
        if not os.path.exists(file_name):
            print("")
            print_info('sleep_x_times_waiting_for_file_presence_on_the_file_system','the file %s does not yet exist on the file system' %file_name )
            count = 0

            while not os.path.exists(file_name):
                if str(sleep_count) != 'wait_til_done':
                    if count < int(sleep_count):
                        print_info('sleep_x_times_waiting_for_file_presence_on_the_file_system','sleeping for 2 seconds')
                        time.sleep(2)
                        count = count + 1
                    else:
                        print_error('leep_x_times_waiting_for_file_presence_on_the_file_system','unable to proceed-- exiting')
                else:
                    time.sleep(2)
                    count = count + 1

            print_info('sleep_x_times_waiting_for_file_presence_on_the_file_system','the file %s now exists on the file system' %file_name )
            print("")

    else:
        print_one_of_the_input_variables_not_set_error_and_exit('sleep_x_times_waiting_for_file_presence_on_the_file_system',file_name,sleep_count)


def convert_from_uncode_to_utf_8(input):
    if isinstance(input, dict):
        return {convert_from_uncode_to_utf_8(key): convert_from_uncode_to_utf_8(value) for key, value in list(input.items())}
    elif isinstance(input, list):
        return [convert_from_uncode_to_utf_8(element) for element in input]
    elif isinstance(input, str):
        return input.encode('utf-8')
    else:
        return input


def utilities_decode_list(data):
    rv = []

    if data:
        if len(data) > 0:
            for item in data:
                if isinstance(item, str):
                    item = item.encode('utf-8')
                elif isinstance(item, list):
                    item = utilities_decode_list(item)
                elif isinstance(item, dict):
                    item = utilities_decode_dict(item)
                rv.append(item)

        elif len(data) == 0:
            print("")
            print_info('utilities_decode_list','input list is empty')
            print("")
    else:
        print_input_variable_not_set_error_and_exit('utilities_decode_list','data')

    return rv

def decode_list(data):
    rv = []

    if data:    
        if len(data) > 0:
            for item in data:
                if isinstance(item, str):
                    item = item.encode('utf-8')
                elif isinstance(item, list):
                    item = utilities_decode_list(item)
                elif isinstance(item, dict):
                    item = utilities_decode_dict(item)

                rv.append(item)

        elif len(data) == 0:
            print("")
            print_info('decode_list','ERROR, input list is empty')
            print("")

    else:
        print_input_variable_not_set_error_and_exit('decode_list','data')

    return rv

def get_key_from_dictionary_given_value(some_dictionary,some_value):

    if some_dictionary and some_value:
        for key, value in list(some_dictionary.items()):
            if some_value == value:
             return key

    else:
        print_one_of_the_input_variables_not_set_error_and_exit('get_key_from_dictionary_given_value',some_dictionary,some_value)

    return None

def utilities_decode_dict(data):
    rv = {}

    if data:
        if len(data) > 0:
            for key, value in list(data.items()):
                if isinstance(key, str):
                    key = key.encode('utf-8')
                if isinstance(value, str):
                    value = value.encode('utf-8')
                elif isinstance(value, list):
                    value = utilities_decode_list(value)
                elif isinstance(value, dict):
                    value = utilities_decode_dict(value)
                rv[key] = value

        else:
            print("")
            print_info('utilities_decode_dict','ERROR: empty dictionary"')
            print("")

    else:
        print_input_variable_not_set_error_and_exit('utilities_decode_dict','data')

    return rv

def decode_dict(data):
    rv = {}

    if data:
        if len(data) > 0:
            for key, value in list(data.items()):
                if isinstance(key, str):
                    key = key.encode('utf-8')
                if isinstance(value, str):
                    value = value.encode('utf-8')
                elif isinstance(value, list):
                    value = utilities_decode_list(value)
                elif isinstance(value, dict):
                    value = utilities_decode_dict(value)
                rv[key] = value

        else:
            print("")
            print_info('decode_dict','ERROR: empty dictionary')
            print("")

    else:
        print_input_variable_not_set_error_and_exit('decode_dict','data')

    return rv

def sysdate(format=None):
    if not format:
        format = "%a %b %d, %Y, %H:%M"

    now = datetime.datetime.now()
    today = now.strftime(format)

    return today


def create_portal_from_external_script(some_value):
    dataLoad_path = None

    if some_value:
        if not os.environ.get('DATALOAD_HOME'):
            print_error('create_portal_from_external_script','this python function needed is needed to be source with config/env3.0.csh')
        else:
            dataLoad_path = os.environ.get('DATALOAD_HOME')

        create_portal_tool = str(dataLoad_path) + '/custom_python_modules/docs/create_portal.py'

        cmd = str(create_portal_tool) + " " + str(some_value)
        print("")
        print_info('create_portal_from_external_script','running the bash command %s' %cmd)
        print("")
        run_bash_command(cmd)

    else:
        print_input_variable_not_set_error_and_exit('create_portal_from_external_script','some_value')

    return some_value

def is_there_another_instance_running(script_name):
    # ps -efa | grep <invoking_script's name> | grep python | wc -l
    #
    # on the sript calling this function, 
    # (1) script_name = os.path.basename(__file__)
    # (2) is_there_another_instance_running(script_name)
    #
    # WHEN INVOKING THIS FUNCTION, check the count from the script that
    # calls this script since it might match more than the script_name
    # since the invoking script name also matches it-- example:
    # str(BASE_DIR) + /bin/taxon_oid_assignment.py checks for
    # is_there_another_instance_running w/string "taxon_oid_assignment.py"
    # and it will alway be 1 or more since theinvoking script will also be in theprocess table
    ps, grep1, grep2, wc, wc_count = None, None, None, None, None

    if script_name:
        hostname = socket.gethostname()
        cmd = "ps -efa | grep %s |grep -v grep | wc -l" %script_name
        wc_count = subprocess.check_output(cmd, shell=True, encoding='utf8').strip('\n')
        if wc_count != None:
            if int(wc_count) == 0:
                print("")
                print_info('is_there_another_instance_running','no other copy of %s running on the host %s' %(script_name, hostname))
                print("")
                return wc_count, "sole_script_intances_running on hostname %s" %hostname

            elif int(wc_count) > 0:
                print("")
                print_info('is_there_another_instance_running','%s instances of %s running on the host %s' %(wc_count, script_name, hostname))
                print("")
                return wc_count, "other_script_intances_running"

        else:
            print_error('is_there_another_instance_running','this is an issue, the wc_count is %s' %wc_count)
    else:
        print_input_variable_not_set_error_and_exit('is_there_another_instance_running','script_name')

    return wc_count, "do not know if other copies of %s are running on this host" %script_name

def run_cypher_on_string(some_string, seed, action):
    import cryptocode

    transformed_string = None

    if some_string and seed and action:
        if action == 'decrypt' or action == 'encrypt':
            if action == 'encrypt':
                transformed_string = cryptocode.encrypt(some_string,seed)
            elif action == 'decrypt':
                transformed_string = cryptocode.decrypt(some_string, seed)
        else:
            print_error('run_cypher_on_string', 'the only action is either encrypt or decrypt, not %s' %action)
    else:
        print_one_of_the_input_variables_not_set_error_and_exit('run_cypher_on_string',some_string, seed, action)

    return transformed_string

def get_nersc_passwd_opt_combo(user_name):
    # for a given user_name, it will provide the passwd + OPT combo
    # avoid calling this function more than 5 min per invocation.

    import pyotp

    key_file = "/home/pphajek__lbl.gov/.ssh/info"
    password_combo, decrypt_key, secret_enc, base32secret, passwd_enc, passwd = None, None, None, None, None, None

    if user_name:
        if os.path.exists(key_file):
            key_dict = load_dictionary_from_file(key_file, '\t')
            #print_dict(key_dict,'key_dict')
            if user_name == 'phajek':
                if key_dict.get('A:'):
                    secret_enc = key_dict[ 'A:' ]
                else:
                    print_not_set_error_and_exit('get_nersc_passwd_opt_combo','A:')
                if key_dict.get('C:'):
                    passwd_enc = key_dict[ 'C:' ]
                else:
                    print_not_set_error_and_exit('get_nersc_passwd_opt_combo','C:')

            if key_dict.get('B:'):
                decrypt_key = key_dict[ 'B:' ]
            else:
                print_not_set_error_and_exit('get_nersc_passwd_opt_combo','B:')

            if secret_enc and passwd_enc and decrypt_key:
                passwd = run_cypher_on_string(passwd_enc, str(decrypt_key), 'decrypt')
                base32secret = run_cypher_on_string(secret_enc, str(decrypt_key), 'decrypt')
                otp_uri = pyotp.totp.TOTP(base32secret).provisioning_uri( "pphajek@lbl.gov", issuer_name="LinOTP")
                otp = pyotp.TOTP(base32secret)
                otp_key = otp.now()
                password_combo = str(passwd) + str(otp_key)
            else:
                print_one_of_the_variables_not_set_error_and_exit('get_nersc_passwd_opt_combo',secret_enc,passwd_enc,decrypt_key)
        else:
            print_path_does_not_exist_error_and_exit('get_nersc_passwd_opt_combo',key_file)
    else:
        print_input_variable_not_set_error_and_exit('get_nersc_passwd_opt_combo','user_name')

    return password_combo

def rsync_from_host_using_file_listing(user,passwd,src_host,dest_base_dir,file_list_name,timeout=3600):
    # rsync will read a file containing a list of remote file paths (one file path per line) located at remote 
    # host src_host and rsync the files to the host calling this function, i.e.
    #
    # we want to get the files from  
    # src_host:/global/cfs/cdirs/img/annotated_submissions/261760/Ga0506620_functional_annotation.gff
    # src_host:/global/cfs/cdirs/img/annotated_submissions/261760/Ga0506620_gene_phylogeny.tsv
    # src_host:/global/cfs/cdirs/img/annotated_submissions/261760/Ga0506620_ko_ec.gff
    # 
    # the file_list_name path should contain the 3 lines:
    # /global/cfs/cdirs/img/annotated_submissions/261760/Ga0506620_functional_annotation.gff
    # /global/cfs/cdirs/img/annotated_submissions/261760/Ga0506620_gene_phylogeny.tsv
    # /global/cfs/cdirs/img/annotated_submissions/261760/Ga0506620_ko_ec.gff
    # 
    # and after this function is called, the files will be located at:
    # str(dest_base_dir) + '/global/cfs/cdirs/img/annotated_submissions/261760/Ga0506620_functional_annotation.gff'
    # str(dest_base_dir) + '/global/cfs/cdirs/img/annotated_submissions/261760/Ga0506620_gene_phylogeny.tsv'
    # str(dest_base_dir) + '/global/cfs/cdirs/img/annotated_submissions/261760/Ga0506620_ko_ec.gff'

    coverage_file = None

    if user and passwd and src_host and dest_base_dir and file_list_name:
        if os.path.exists(file_list_name):
            #file_list_name = str(dataLoad_path) + '/tmp/scratch_dir/file_list.txt'
            rsync_cmd = '/usr/bin/sshpass -p \'%s\' /usr/bin/rsync -av -e \'ssh -q -o StrictHostKeyChecking=no\' --files-from=%s %s@%s:/ %s' %(passwd,file_list_name,user,src_host,dest_base_dir)
            print(rsync_cmd)
            os.system(rsync_cmd)
            #line = '/usr/bin/sshpass -p \'%s\' /usr/bin/rsync -av -e \'ssh -q -o StrictHostKeyChecking=no\' %s@%s:%s %s/' %(passwd,user,src_host,coverage_file,dst_dir)
        else:
            print_path_does_not_exist_error_and_exit('rsync_from_host_using_file_listing',file_list_name)
    else:
        print_one_of_the_input_variables_not_set_error_and_exit('rsync_from_host_using_file_listing',user,passwd,src_host,dest_base_dir,file_list_name)

def rsync_processed_merfs_content_to_nersc(taxon_oid):
    # first generate a list of files and directories to be rsync'ed to NERSC'd dna volume
    # then rsync content
    
    if taxon_oid:
        if os.environ.get('DATALOAD_HOME'):
            dataLoad_path = os.environ['DATALOAD_HOME']
            file_list_name = str(dataLoad_path) + '/tmp/scratch_dir/file_list_to_nersc.txt'
            if os.path.exists(file_list_name):
                os.remove(file_list_name)

            src_merfs_dir = '/projectdirs/img_web_data_merfs/' + str(taxon_oid) + '/'
            src_taxon_tarball = '/projectdirs/img_web_data_merfs/download/' + str(taxon_oid) + '.tar.gz'
            src_phylodist_dir = '/projectdirs/img_web/img_web_data_ava/phyloDist/' + str(taxon_oid) + '/'
            content_list = [ src_merfs_dir, src_taxon_tarball, src_phylodist_dir ]
            if os.path.exists(src_merfs_dir) + os.path.exists(src_taxon_tarball) + os.path.exists(src_phylodist_dir):
                dest_merfs_dir = '/global/dna/projectdirs/microbial/img_web_data_merfs/' + str(taxon_oid)
                dest_phylodist_dir = '/global/dna/projectdirs/microbial/img_web_data_merfs/phyloDist/' + str(taxon_oid)
                dest_taxon_tarball = '/global/dna/projectdirs/microbial/img_web_data_ava/download/'+ str(taxon_oid)

                with open(file_list_name, 'w') as f:
                    for file_path in content_list:
                        line = str(file_path) + '\n'
                        f.write(str(line))

                if os.path.exists(file_list_name):
                    user = 'phajek'
                    dest_host = 'dtn04.nersc.gov'
                    dest_dir = '/global/dna/projectdirs/microbial/img_web_data/mer_submisisons_unsorted'
                    passwd = get_nersc_passwd_opt_combo(user)
                    if passwd:
                        rsync_cmd = '/usr/bin/sshpass -p \'%s\' /usr/bin/rsync -avr -e \'ssh -q -o StrictHostKeyChecking=no\' --files-from=%s / %s@%s:%s' %(passwd,file_list_name,user,dest_host,dest_dir)
                        print(rsync_cmd)
                        os.system(rsync_cmd)
                    else:
                        print_not_set_error_and_exit('rsync_processed_merfs_content_to_nersc','rsync_cmd')
                else:
                    print_path_does_not_exist_error_and_exit('rsync_content_to_nersc.py', file_list_name)
            else:
                print_path_does_not_exist_error_and_exit('rsync_content_to_nersc.py', src_merfs_dir)
                print_path_does_not_exist_error_and_exit('rsync_content_to_nersc.py', src_taxon_tarball)
                print_path_does_not_exist_error_and_exit('rsync_content_to_nersc.py', src_phylodist_dir)
        else:
            print_not_set_error_and_exit('rsync_processed_merfs_content_to_nersc','DATALOAD_HOME')
    else:
        print_input_variable_not_set_error_and_exit('rsync_processed_merfs_content_to_nersc','taxon_oid')
       
def get_igb_password(user_name):
    passwd_enc, seed, password = None, None, None

    if user_name:
        key_file = '/global/homes/i/img/.ssh/info'
        if os.path.exists(key_file):
            key_dict = load_dictionary_from_file(key_file, '\t')
            if key_dict:
                if user_name == 'phajek':
                    if key_dict.get('A:'):
                        passwd_enc = str(key_dict[ 'A:' ])
                    else:
                        print_not_set_error_and_exit('get_igb_password','A:')

                    if key_dict.get('B:'):
                        seed = str(key_dict[ 'B:' ])
                    else:
                        print_not_set_error_and_exit('get_igb_password','B:')
                else:
                    print_error('get_igb_password','so far, only user phajek has be configured to retrieve the password')

                password = run_cypher_on_string(passwd_enc, seed, 'decrypt')
            else:
                print_not_set_error_and_exit('get_igb_password','key_dict')
        else:
            print_path_does_not_exist_error_and_exit('get_igb_password',key_file)
    else:
        print_input_variable_not_set_error_and_exit('get_igb_password','user_name')
        
    return password

def main(base_dir, filename):
    submission_id, taxon_oid, full_file_path, taxon_oid_dir = None, None, None, None
    img, metadata = {}, {}

    userName = 'pphajek'
    full_file_path = str(base_dir) + '/' + str(filename)

    match = re.match('(.*)_.*\.tar\.gz', filename)
    if match:
            taxon_oid = match.group(1)
            dsn_tns = cx_Oracle.makedsn('oracle-vm-1.jgi.lbl.gov', '1521', service_name='IMGIPRDPDB')
            connection = cx_Oracle.connect(user=r'imgsg_dev', password='Tuesday', dsn=dsn_tns)
            cursor = connection.cursor()

            sql = "select submission_id from submission where img_taxon_oid=%s" %taxon_oid
            print(("Finding submission_id %s" %sql))
            cursor.execute(sql)
            for row in cursor.fetchall():
                if row:
                    submission_id = int(row[0])
                    print(("The submission_id: %s" %submission_id))
                    img[ 'taxon_oid' ] = int(taxon_oid)
                    img[ 'submission id' ] =int(submission_id)
                    img[ 'full file path' ] = full_file_path
                    gold_data,img,metadata = get_submission_metadata_data(img,metadata)
                    portal_id, portal_status = find_portal_id(img,metadata)
                    submit_to_jamo(gold_data,img,metadata)
                else:
                    print(("**** No submission id associated with taxon oid %s" %taxon_oid))


if __name__ == "__main__":
    if int(len(sys.argv))== 2:
        dataLoad_path, base_dir, filename = None, None, None
    else:
        sys.exit(1)
