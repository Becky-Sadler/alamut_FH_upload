import pandas as pd
import numpy as np 

# create a dataframe from the data (after converting it to a csv file in excel)
df = pd.read_csv('test.csv')
df = df.dropna(how = 'all')

# create empty dataframe for alamut ready variants
cNomen_columns = ["Gene", "Transcript", "Variant", "Pathogenic", "Patient ID", "Family ID", "Phenotype", "Comment"]
cNomen_df = pd.DataFrame(columns = cNomen_columns)


# Loop to process the data into two different tab delimited .txt files for upload to Alamut. 
for index, row in df.iterrows():
        if row['GOSH ID'] == 'FH Email Group Variants' or row['GOSH ID'] == 'Referring Lab':
                continue
        else:
                Transcript, cNomen = row['acmg_cdna'].split(':')
                PatientID = row['GOSH ID']
                Gene = row['gene']
                if row['acmg_classification'] == '1' or row['acmg_classification'] == '2' or row['acmg_classification'] == '3' or row['acmg_classification'] == '4' or row['acmg_classification'] == '5':
                        Pathogenic = 'Class ' + str(row['acmg_classification'])
                else:
                        Pathogenic = 'Class 3'
                # removing new lines from the comment column
                string = str(row['comment']).replace('\n', ' ').replace('\r', '')
                # Combining evidence column with the comment column to add to the general Alamut comment column (separated with a : as cannot do a new line separation)
                Comment = str(row['criteria']) + ' : ' + string + ' : Seen :' + str(row['seen'])  
                # Append this record to the cNomen dataframe 
                cNomen_df = cNomen_df.append({'Gene': Gene, 'Transcript': Transcript, 'Variant' : cNomen, 'Pathogenic': Pathogenic, 'Patient ID' : PatientID, 'Comment' : Comment}, ignore_index=True)

APOB_df = cNomen_df.loc[(cNomen_df['Gene'] == 'APOB')]
APOB_df.to_csv('FH_APOB.txt', sep='\t', index = False)

APOE_df = cNomen_df.loc[(cNomen_df['Gene'] == 'APOE')]
APOE_df.to_csv('FH_APOE.txt', sep='\t', index = False)

LDLR_df = cNomen_df.loc[(cNomen_df['Gene'] == 'LDLR')]
LDLR_df.to_csv('FH_LDLR.txt', sep='\t', index = False)

LDLRAP1_df = cNomen_df.loc[(cNomen_df['Gene'] == 'LDLRAP1')]
LDLRAP1_df.to_csv('FH_LDLRAP1.txt', sep='\t', index = False)

PCSK9_df = cNomen_df.loc[(cNomen_df['Gene'] == 'PCSK9')]
PCSK9_df.to_csv('FH_PCSK9.txt', sep='\t', index = False)
