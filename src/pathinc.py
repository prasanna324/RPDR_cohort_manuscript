
def load_RPDR_path_multiple(dir_data, fname, delimiter='|', datetime_col='Report_Date_Time'):
    ''' load_RPDR_path_multiple(dir_data, fname, delimiter='\t', datetime_col='Report_Date_Time'):
        Sequentially loads all files from RPDR data dump when output is split. 
        
        1. Starts in dir_data (should have trailing slash), grabs all sub-folders' names automatically, then sequentially loads: dir_data/[sub-folders]/fname (where fname is the name of the file)
        * Note for whatever reason, on a multiple-split file dump from RPDR the labs, demographics, etc files are all named the exact same, just in different zips
        2. Calls the traditional load path function on each file
        3. Concatenates all results and returns 1 DF
        
        See load_native_data for remainder of parameters which are passed to that function
        
        '''
    import os
    import pandas as pd
    
    # get list of subdirectories
    subdirectories = [x[0] for x in os.walk(dir_data)][1:]
    
    first=True
    # for each subdir, use the traditional load function to load data and concat
    for subdir in subdirectories:
        path_to_path=subdir+'/'+fname
        path_to_path_processed = path_to_path.replace('.txt','_multiline_corrected.txt')
        if (not os.path.exists(path_to_path)) & (not os.path.exists(path_to_path_processed)):
            continue
        path = load_RPDR_path(path=path_to_path,
                              delimiter=delimiter,
                              datetime_col=datetime_col)
        
        if first==True:
            concat_pd = path
            first=False
        else:
            concat_pd=pd.concat([concat_pd, path],ignore_index=True)
    
    return concat_pd

def load_RPDR_path(path,delimiter='|', datetime_col='Report_Date_Time'):
    ''' load_RPDR_path(path,string_format=False,delimiter='|', datetime_col='Report_Date_Time')
    DESC: loads an RPDR pathology file to pandas 
    1. removes deleted, canceled and in process path reports
    2. removes duplicate path reports, keeping the first entry
    3. resets index to 'Report_Number' column
    4. converts datetime_col to pd.DateTime format
    
    PARAMETERS:
    path: path to csv file or other text delimited file
    delimiter: delimiter for path file
    datetime_col: column name containing date/time information for each path report
    
    returns: pandas dataframe containing path information
    
    WARNINGS:
    1. Current function automatically searches for path + 'multiline_corrected', *if present* it assumes that is the correct 
        file. E.g., path='/data/path.txt', it searches for '/data/path_multiline_corrected.txt'.
    2. It will not overwrite this file if present
    
    TO-DO:
    1. Update path report selection when there are near-duplicates. These are almost always of the form 'Final' and 'Updated', 
        of which, when duplicated, we should take the 'Updated' report. I have looked at about 10 examples so far (Marc)
    '''
    import pandas as pd
    import os.path
    from os import path as os_path
    
    ## default Report_Text format is multi-line. If string_format==True, convert multi-line text by double-quoting all quotes
    ##  (which allows read_csv to read through by default), and enclosing full report in single-quotes
    write_path = path.replace('.txt','_multiline_corrected.txt')
    if os_path.exists(write_path)==False:
        print('Reformatting path file to allow multi-line report text to be readable, saving as : {}'.format(write_path))
        f = open(path,'r')
        filedata = f.read()
        f.close()
        newdata = filedata.replace('"', '""')
        ## Mod 1 ##
        # Replacing Accession number with PAT to deal with Report_Text's that doesn't contain the word "Accession number"
        newdata = newdata.replace('|PAT|', '|PAT|"')
        ## Mod 1 ##
        newdata = newdata.replace('[report_end]', '[report_end]"')
        f2 = open(write_path,'w')
        f2.write(newdata)
        f2.close()
        
    path = write_path
    
    print('Reading from : ' + path)
    path_df = pd.read_csv(path, sep=delimiter, dtype=str)
    
    # unique_identifier. Requires some explanation.. report_number SHOULD be unique (it is literally the accession
    #  number for finding a block of tissue), however it is NOT in some specific examples of NSMC and MGH. E.g., 
    #  in all_RPDR_path, report numbers  : 
    # 'S11-5510', 'S14-3070', 'S17-25856', 'S13-2400', 'S14-9841', 'S13-3071',
    #    'S12-6506', 'S14-218', 'S12-3414', 'S13-32', 'S13-6114', 'S12-2481',
    #    'S13-1077', 'S14-41', 'S13-8207', 'S12-10612', 'S12-5964', 'S10-9798',
    #    'S09-10295', 'S16-15842', 'S14-8374', 'S12-3350', 'S12-6495', 'S14-785',
    #    'S16-183'
    #  all give duplicates for NSMC & MGH. And they're *different* patients, totally diff reports
    #  solution: append EMPI (unique per patient) to report_number to screen out these cases, operate on this concatenated report id
    path_df['unique_report_id'] = path_df.apply(lambda x: str(x.EMPI) + '_' + str(x.Report_Number),axis=1)
    
    # remove deleted, cancelled, and in-process path reports
    # bad_reports = (path_df.Report_Status == 'Deleted') | (path_df.Report_Status == 'Cancelled') | (path_df.Report_Status == 'In Process') | (path_df.Report_Status == 'Preliminary') | (path_df.Report_Status == 'Hold')  | (path_df.Report_Status == 'Pending')
    bad_statuses=['Deleted', 'Cancelled', 'In Process', 'Preliminary', 'Hold', 'Pending','Unknown','In Revision', 'Not Verified']
    bad_reports = path_df.Report_Status.isin(bad_statuses)
    
    # drop rows with any of these features
    path_df2 = path_df[~bad_reports].copy()
    
    ## Mod 2 ##
    # Drop duplicates for 'unique_report_id' cases based on length of Result_Text instead of dropping the first observed case
    path_df2['report_len'] = path_df2['Report_Text'].str.len()
    
    path_df2 = (path_df2
      .sort_values(['unique_report_id', 'report_len'])
      .drop_duplicates(subset=['unique_report_id'], keep='last', inplace=False, ignore_index=False)
      .drop(columns=['report_len'])
      .sort_index()
     )
    ## Mod 2 ##
    
    # now set index to unique_report_id
    path_df2.set_index(keys='unique_report_id', inplace=True, verify_integrity=True)
    
    # set date time as datetime
    path_df2['datetime'] = pd.to_datetime(path_df2.loc[:,datetime_col])

    return path_df2

def truncate_finaldx(pathdf, update=True):
    ''' truncate_finaldx(pathdf, update=True)
    DESC: For path reports, find the 'final diagnosis' line of the path report, remove everything preceding. Parse 
     to extract whether or not there is a final diagnosis line, what that line is, and the full report text (following this line)
    
    PARAMETERS:
    pathdf: pathology dataframe from load_RPDR_path
    update: return preprocessed path as a new df or update pathdf
    
    RETURNS: pathdf or new path dataframe with only MGH, BWH path reports (depending on update bool) with three new columns:
    ['has_final_diagnosis'] = did the function find a line of text that it thinks contains final diagnosis line?
    ['final_diagnosis_line'] = if has_final_diagnosis == True, then what is the final diagnosis line?
    ['Report_Text'] = report text after removing everything above final diagnosis
    '''
    import re
    
    # truncate to only final diagnosis
    
    # first get only MGH values
    #fil_mgh = pathdf.MRN_Type == 'MGH'
    df_path = pathdf.copy()
    
    # FINAL DIAGNOSIS LINE FINDER
    num_reports = df_path.shape[0] # num rows of df_path
    has_final_diagnosis_col = []
    final_diagnosis_line_col = []
    trunc_path_col = []
    print('Truncating to only final diagnosis...')
    for i in range(0,num_reports):
        # extract path report for this entry
        report_text = df_path.iloc[i,:].Report_Text
        site = df_path.iloc[i,:].MRN_Type
        description = df_path.iloc[i,:].Report_Description
        # split by newline character
        text_by_line = report_text.split('\n')

        has_final_diagnosis = False
        final_diagnosis_line = ''
        trunc_path_text = report_text
        non_excl = not bool(re.search(r'\bfetus\b|\bfetopsy\b|\bautopsy\b', report_text.lower()))

        # go line-by-line and perform some checks
        j=0
        
        if site in ['MGH','BWH','NWH','FH','NSM'] and non_excl:
            for line in text_by_line:
                lower_line = line.lower().strip()
                # capture situation where a line contains liver, biopsy; note will only grab first instance then short circuit
                
                if (has_final_diagnosis==False) and ((
                            ('Final' in line or 'Pathologic' in line or
                            'FINAL' in line or 'PATHOLOGIC' in line) and 'diagnosis' in lower_line) or
                            (lower_line.startswith('diagnosis:')) or
                            (line.strip()=='Diagnosis') or
                            ('FINAL REPORT' in line) or
                            (line.startswith('SPECIMEN(S):'))) and not (
                                            'amend' in lower_line or
                                            'final diagnosis by' in lower_line or
                                            'clinical data' in lower_line or
                                            'cytogenetic' in lower_line or
                                            'reason for' in lower_line or
                                            'gross description' in lower_line or
                                            'original' in lower_line or
                                            'epithelial' in lower_line or
                                            'unsatisfactory for evaluation' in lower_line or
                                            'clinical history' in lower_line):
                    has_final_diagnosis = True
                    final_diagnosis_line = line
                    trunc_path_text = '\n'.join(text_by_line[j:]) # should be a list of this line and all subsequent lines

                j=j+1
                
        has_final_diagnosis_col.append(has_final_diagnosis)
        final_diagnosis_line_col.append(final_diagnosis_line)
        # either returns the original report or the truncated form if it has a final diagnosis to truncate at
        trunc_path_col.append(trunc_path_text)

    df_path['has_final_diagnosis'] = has_final_diagnosis_col
    df_path['final_diagnosis_line'] = final_diagnosis_line_col
    df_path['Report_Text_Raw'] = df_path['Report_Text']
    df_path['Report_Text'] = trunc_path_col


    if update:
        # re-merge with original data
        print('Updating input path dataframe with truncated path reports')
        pathdf['has_final_diagnosis'] = False
        pathdf['final_diagnosis_line'] = ''
        pathdf['Report_Text_Raw'] = pathdf['Report_Text']
        pathdf.update(df_path)
        return_df = pathdf.copy()
    else:
        # return this mgh, bwh path only file
        print('Returning entries with truncated path reports')
        return_df = df_path
        
        print('Done. | Status: ' + str(df_path[df_path.has_final_diagnosis==True].shape[0]) + ' reports with a final diagnosis, ' 
              + str(df_path[df_path.has_final_diagnosis==False].shape[0]) + ' reports with no final diagnosis')
    
    return return_df


def truncate_lower(pathdf, update=True, only_finaldx=True):
    
    #print('Filtering only MGH, BWH path reports...')
    fil_subset = pathdf.MRN_Type.isin(['MGH', 'BWH','NWH','FH','NSM'])
    df_path = pathdf[fil_subset].copy()
    
    if only_finaldx:
        # check the column exists first:
        if 'has_final_diagnosis' in df_path.columns.tolist():
            fil_finaldx_trunc = df_path.has_final_diagnosis == True
            df_path = df_path[fil_finaldx_trunc]
        else:
            print('The flag *only_finaldx=True* was passed, however truncate_finaldx() has not been called. Aborting...')
            return None

    num_reports = df_path.shape[0]
    has_lowersec_col = []
    lowersec_line_col = []
    lowersec_start_LAFD_col = []
    trunc_path_col = []
    
    has_addendum_col = []
    
    for i in range(0,num_reports):
        # extract path report for this entry
        report_text = df_path.iloc[i,:].Report_Text
        # split by newline character
        text_by_line = report_text.split('\n')
        
        has_lowersec = False
        lowersec_line = ''
        lowersec_start_LAFD = -1
        trunc_path_text = report_text
        
        has_addendum = False
        addendum_start_LAFD = -1
        addendum_end_LAFD = len(text_by_line)
        
        # go line-by-line and perform some checks
        j=0
        for line in text_by_line:
            lower_line = line.lower()
            
            trimmed_line = remove_extra_spaces(line)
            trimmed_lower_line = remove_extra_spaces(lower_line)
            
            if ('CLINICAL DATA' in trimmed_line or
                                        'by his/her signature below' in trimmed_lower_line or
                                        'Final Diagnosis by' in trimmed_line or
                                        'electronically signed out' in trimmed_lower_line or
                                        'diagnosis by:' in trimmed_lower_line or
                                        'gross description:' in trimmed_lower_line or
                                        'o. r. consultation report:' in trimmed_lower_line or
                                        'Reports to:' in trimmed_line or
                                        'SPECIMEN TYPE:' in trimmed_line or
                                        'clinical diagnosis & history' in lower_line or
                                        'pathology dept report date/time' in lower_line or
                                        'Clinical History:'==line or
                                        ('addend' in lower_line and 'by' in lower_line)):
                if (has_addendum==False and has_lowersec==False):
                    has_lowersec = True
                    lowersec_line = line
                    lowersec_start_LAFD = j
                    addendum_start_LAFD = j
                if has_addendum==True:
                    addendum_end_LAFD = j
                
            if has_addendum==False and (
                (('ADDENDUM' in trimmed_line or 'Addendum' in trimmed_line) and (':' in trimmed_line)) or
                ('ADDENDUM'==trimmed_line) or ('Addendum'==trimmed_line)
            ):
               
                has_addendum = True
                addendum_start_LAFD = j
                
            j=j+1
        
        if has_lowersec==True:
            
            if has_addendum==True:
                trunc_path_text = '\n'.join(text_by_line[:lowersec_start_LAFD]+text_by_line[addendum_start_LAFD:addendum_end_LAFD])
            else:
                trunc_path_text = '\n'.join(text_by_line[:lowersec_start_LAFD])
            

        has_lowersec_col.append(has_lowersec)
        lowersec_line_col.append(lowersec_line)
        lowersec_start_LAFD_col.append(lowersec_start_LAFD)
        trunc_path_col.append(trunc_path_text)
        has_addendum_col.append(has_addendum)
        
    df_path['has_lowersec'] = has_lowersec_col
    df_path['lowersec_line'] = lowersec_line_col
    df_path['lowersec_start_LAFD'] = lowersec_start_LAFD_col
    df_path['Report_Text'] = trunc_path_col
    df_path['has_addendum'] = has_addendum_col
    
    
    if update:
        # re-merge with original data
        print('Updating input path dataframe with truncated path reports')
        pathdf['has_lowersec'] = False
        pathdf['lowersec_line'] = ''
        pathdf['lowersec_start_LAFD'] = -1
        pathdf['has_addendum'] = False
        pathdf.update(df_path)
        return_df = pathdf.copy()
    else:
        # return this mgh path only file
        print('Returning entries with truncated path reports')
        return_df = df_path
    
    return return_df


def is_liver_biopsy(pathdf, update=True, only_finaldx=True): 
    ''' is_liver_biopsy(pathdf, update=True, only_finaldx=True)
    DESC: For path reports, determine whether this pathology report is from a liver biopsy. 
    REQUIRES: a column named Report_Text containing the full text of the path report; if only_finaldx is true, must have run
     truncate_finaldx() first 
    
    PARAMETERS:
    pathdf: pathology dataframe from load_RPDR_path (or subsequent modification)
    update: return preprocessed path as a new df or update pathdf results
    only_finaldx: bool, if True, only update rows where has_final_diagnosis==True
    
    RETURNS: pathdf or new path dataframe with preprocessed path reports (depending on update bool) with three new columns:
    ['is_liver_biopsy'] = bool, is this a liver biopsy or not?
    ['is_liver_biopsy_line'] = text of liver biopsy line if above true
    ['liver_biopsy_LAFD'] = LAFD=Lines After Final Diagnosis, ie if final diagnosis line on line 12, and liver biopsy line 16, =16-12
     This is to help screen for oddities, as there is a typical, rough number of lines between these entries.
    '''    
    
    import re
    
    #print('Filtering only MGH, BWH path reports...')
    #fil_mgh = pathdf.MRN_Type == 'MGH'
    df_path = pathdf.copy()
    
    if only_finaldx:
        # check the column exists first:
        if 'has_final_diagnosis' in df_path.columns.tolist():
            fil_finaldx_trunc = df_path.has_final_diagnosis == True
            df_path = df_path[fil_finaldx_trunc]
        else:
            print('The flag *only_finaldx=True* was passed, however truncate_finaldx() has not been called. Aborting...')
            return None
    
    num_reports = df_path.shape[0]
    is_liver_biopsy_col = []
    is_liver_biopsy_line_col = []
    liver_biopsy_LAFD_col = []
    for i in range(0,num_reports):
        # extract path report for this entry
        report_text = df_path.iloc[i,:].Report_Text
        # split by newline character
        text_by_line = report_text.split('\n')

        is_liver_biopsy = False
        is_liver_biopsy_line = ''
        liver_biopsy_LAFD = -1

        # go line-by-line and perform some checks
        j=0
        for line in text_by_line:
            lower_line = line.lower().lstrip()
            # capture situation where a line contains liver, biopsy; note will only grab first instance then short circuit
            if is_liver_biopsy==False and ((bool(re.search(r'\bliver\b|\bhepatic\b',lower_line)) and (
                        bool(re.search(r"""\bbiopsy\b|\bbiopsies\b|\blobe\b|\brandom\b|
                        |\bnon-focal\b|\bnonfocal\b|\bnon focal\b|\bcore\b|bcores\b|
                        |\bspecimen\b|\bmass\b|\blesion\b|\btissue\b|\bparenchyma\b|
                        |\boperation\b|\bsegment\b|\bsegments\b|\bexcision\b|\bnodule\b|
                        |\bexplant\b|\bresection\b|\bhepatectomy\b|
                        |\blobectomy\b|\bnative\b""",lower_line)))) or 
                               ('hepatectomy' in lower_line and 
                                bool(re.search(r'\bleft\b|\bright\b|\bpartial\b|\bcentral\b',lower_line))) or (
                                        bool(re.search(r'^\s*[a-f][.].*\bliver\b',lower_line)) or
                                        bool(re.search(r'^\s*[a-f]/1.*\bliver\b',lower_line)) or
                                        bool(re.search(r'^\s*[a-f][.|)].*\bliver\b.*[(.*)]:',lower_line)) or
                                        bool(re.search(r'.*\bliver\b.*[:]', lower_line)))) and not (
                                                'colon' in lower_line or 
                                                'renal' in lower_line or
                                                'kidney' in lower_line or
                                                'flexure' in lower_line or
                                                'flex' in lower_line or
                                                'metastatic' in lower_line or
                                                'glomerulosclerosis' in lower_line or
                                                'lymph' in lower_line or
                                                'brush' in lower_line or
                                                'artery' in lower_line or
                                                'fluid' in lower_line or
                                                'cautery artifact' in lower_line or
                                                'duct' in lower_line or
                                                'clinical data' in lower_line or
                                                'gross description' in lower_line or
                                                'original pathologic' in lower_line or
                                                'fna' in lower_line or
                                                'aspiration' in lower_line or
                                                'bile' in lower_line or
                                                'history' in lower_line or
                                                 lower_line.startswith('specimen:')):
                                            #bool('liver:'==lower_line.strip())
                
                total_lc = sum([1 if x.islower() else 0 for x in line.split()])
                total_uc = sum([1 if x.isupper() else 0 for x in line.split()])
                if not ((total_lc>3) & (total_uc<3)) and total_lc<=8:
                    is_liver_biopsy = True
                    is_liver_biopsy_line = line
                    liver_biopsy_LAFD = j
                
            j=j+1

        is_liver_biopsy_col.append(is_liver_biopsy)
        is_liver_biopsy_line_col.append(is_liver_biopsy_line)
        liver_biopsy_LAFD_col.append(liver_biopsy_LAFD)

    df_path['is_liver_biopsy'] = is_liver_biopsy_col
    df_path['is_liver_biopsy_line'] = is_liver_biopsy_line_col
    df_path['liver_biopsy_LAFD'] = liver_biopsy_LAFD_col

    if update:
        # re-merge with original data
        print('Updating input path dataframe with truncated MGH, BWH path reports')
        pathdf['is_liver_biopsy'] = False
        pathdf['is_liver_biopsy_line'] = ''
        pathdf['liver_biopsy_LAFD'] = -1
        pathdf.update(df_path)
        return_df = pathdf.copy()
    else:
        # return this mgh, bwh path only file
        print('Returning entries with truncated path reports')
        return_df = df_path
        
    print('Done. | Status: ' + str(df_path[df_path.is_liver_biopsy==True].shape[0]) + ' reports with likely liver biopsy, ' 
          + str(df_path[df_path.is_liver_biopsy==False].shape[0]) + ' reports likely not a liver biopsy')
    
    return return_df


def mgh_find_note_start(pathdf, update=True, only_finaldx=True):
    
    if not 'is_liver_biopsy' in pathdf.columns.tolist():
        print('Function requires running is_liver_biopsy() first (uses relative positioning of the biopsy line to call note start)')
        return None
    
    print('Filtering only MGH path reports...')
    fil_mgh = pathdf.MRN_Type == 'MGH'
    mgh_path = pathdf[fil_mgh].copy()
    
    if only_finaldx:
        # check the column exists first:
        if 'has_final_diagnosis' in mgh_path.columns.tolist():
            fil_finaldx_trunc = mgh_path.has_final_diagnosis == True
            mgh_path = mgh_path[fil_finaldx_trunc]
        else:
            print('The flag *only_finaldx=True* was passed, however mgh_truncate_finaldx() has not been called. Aborting...')
            return None
    
    num_reports = mgh_path.shape[0]
    has_note_start_col = []
    note_line_col = []
    note_start_LAFD_col = []
    for i in range(0,num_reports):
        # extract path report for this entry
        report_text = mgh_path.iloc[i,:].Report_Text
        # split by newline character
        text_by_line = report_text.split('\n')

        has_note = False
        has_note_after_biopsy = False
        note_line = ''
        note_start_LAFD = -1

        # go line-by-line and perform some checks
        j=0
        for line in text_by_line:
            lower_line = line.lower()
            # prep for finding the first non-space word in a line
            word_list = remove_extra_spaces(lower_line, return_as_list=True)
            # first conditional: make sure the word list isn't empty
            # second conditional: if we've already found note, don't look further
            # third and beyond conditionals: that the first or second word contains note is best signature of the start of the note section
            # the parentheses are to help it ignore catching note on a new line where it says '(see note)' or (see comment)
            if len(word_list) > 0 and has_note_after_biopsy==False and ((('note' in word_list[0] or 'comment' in word_list[0]) and not ')' in word_list[0]) or 
                                                           (len(word_list) > 1 and (('note' in word_list[1] or 'comment' in word_list[1]) and not ')' in word_list[1]))):
                
                
                # unlike other loop functions above, stipulate the 'best' note comes after the line declaring this to be a liver biopsy
                #  but if no such line materializes, still go with the FIRST line it was seen
                if j > mgh_path.iloc[i,:].liver_biopsy_LAFD:
                    # stop searching, this is a good candidate
                    has_note_after_biopsy = True
                    has_note = True
                    note_line = line
                    note_start_LAFD = j
                # this is the first encounter with a probable note line (has_note is False); keep it unless the if above is triggered with a better one
                elif has_note == False:
                    has_note = True
                    note_line = line
                    note_start_LAFD = j
                
            j=j+1

        has_note_start_col.append(has_note)
        note_line_col.append(note_line)
        note_start_LAFD_col.append(note_start_LAFD)

    mgh_path['has_note_start'] = has_note_start_col
    mgh_path['note_line'] = note_line_col
    mgh_path['note_start_LAFD'] = note_start_LAFD_col

    if update:
        # re-merge with original data
        print('Updating input path dataframe note starts and locations')
        pathdf['has_note_start'] = False
        pathdf['note_line'] = ''
        pathdf['note_start_LAFD'] = -1
        pathdf.update(mgh_path)
        return_df = pathdf.copy()
    else:
        # return this mgh path only file
        print('Returning MGH only entries annotated note starts')
        return_df = mgh_path
        
    print('Done. | Status: ' + str(mgh_path[mgh_path.has_note_start==True].shape[0]) + ' reports with an identifiable note entry, ' 
          + str(mgh_path[mgh_path.has_note_start==False].shape[0]) + ' reports without a note entry')
    
    return return_df

def mgh_extract_finaldx(pathdf, update=True):
    
    if not 'is_liver_biopsy' in pathdf.columns.tolist() and 'has_note_start' in pathdf.columns.tolist():
        print('Function requires running mgh_is_biopsy() and mgh_find_note_start first (uses relative positioning of the biopsy line to call note start)')
        return None
    
    print('Filtering only MGH path reports that are notated as liver biopsies...')
    fil_mgh = pathdf.MRN_Type == 'MGH'
    fil_liverbx = pathdf.is_liver_biopsy == True
    mgh_path = pathdf[fil_mgh & fil_liverbx].copy()
    
    num_reports = mgh_path.shape[0]
    has_finaldx_col = []
    finaldx_text_col = []
    finaldx_LAFD_col = []
    
    for i in range(0,num_reports):
        # extract path report for this entry
        report_text = mgh_path.iloc[i,:].Report_Text
        # split by newline character
        text_by_line = report_text.split('\n')

        has_finaldx = False
        finaldx_text = ''
        finaldx_LAFD = -1

        # grab markers from where biopsy line is noted and where the 'note' line is
        liverbx_LAFD = int(mgh_path.iloc[i,:].liver_biopsy_LAFD) # every entry should have a positive liver biopsy LAFD
        note_LAFD = int(mgh_path.iloc[i,:].note_start_LAFD) # NOT every entry will have a note; that's OK

        # if the order is liver biopsy line >> note, then the final diagnosis usually falls right in between
        if note_LAFD > liverbx_LAFD:
            finaldx_text = ' '.join(text_by_line[liverbx_LAFD+1:note_LAFD])
            # get rid of extra spaces
            finaldx_text = remove_extra_spaces(finaldx_text)
            
            has_finaldx = True
            finaldx_LAFD = liverbx_LAFD+1
        
        has_finaldx_col.append(has_finaldx)
        finaldx_text_col.append(finaldx_text)
        finaldx_LAFD_col.append(finaldx_LAFD)

    mgh_path['has_finaldx_text'] = has_finaldx_col
    mgh_path['finaldx_text'] = finaldx_text_col
    mgh_path['finaldx_LAFD'] = finaldx_LAFD_col

    if update:
        # re-merge with original data
        print('Updating input path dataframe with final diagnosis (short) text')
        pathdf['has_finaldx_text'] = False
        pathdf['finaldx_text'] = ''
        pathdf['finaldx_LAFD'] = -1
        pathdf.update(mgh_path)
        return_df = pathdf.copy()
    else:
        # return this mgh path only file
        print('Returning MGH only entries passing all filters above and annotated with final diagnosis (short) text')
        return_df = mgh_path
        
    print('Done. | Status: ' + str(mgh_path[mgh_path.has_finaldx_text==True].shape[0]) + ' reports with a probable final diagnosis, ' 
          + str(mgh_path[mgh_path.has_finaldx_text==False].shape[0]) + ' reports without a probable final diagnosis')
    
    return return_df

def remove_extra_spaces(input_as_str, return_as_list=False):
    word_list = input_as_str.split(' ')
    while '' in word_list: word_list.remove('')
    
    if return_as_list:
        out = word_list
    else:
        out = ' '.join(word_list)
    
    return out
