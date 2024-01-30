# A main module for FastaHandler created by Hyungtaek Jung

#!/usr/bin/env python3

import sys
import subprocess

usage = '''FastaHandler: Fasta File Manipulation Toolkit
version 1.0.1

Usage: python3 fastahandler.py <module> <parameters>

Modules:
Multi2Single\t| m2s\tConvert a multi-fasta (multiline) into a single-line fasta.
RenameId\t| rid\tRename prefix IDs and headers.
PrefixRename\t| prn\tRename prefix IDs and headers with a user’s input.
PrefixSelectRename\t| psr\tRename prefix IDs and headers with a user’s input (Only).
IdExtract\t| idx\tExtract matched IDs and their corresponding sequences.
IdExtractLocation\t| iel\tExtract matched IDs, locations and their corresponding sequences.
IdExtractLocationMultiple\t| iem\tExtract matched IDs, locations and their corresponding sequences (Multiple).
ReverseComplement\t| rcp\tMake a reverse complement sequence.
FindCountDuplication\t| fcd\tFind and count the duplicated IDs and sequences.
RemoveDuplication\t| rvd\tRemove the duplicated IDs and sequences.
SubsetFasta\t| ssf\tMake a subset of data with a sequence length filter.
ExtractPattern\t| xpt\tMake a subset of data with find, filter and extract.
EachFastaStats\t| efs\tGenerate each line fasta statistic for a multi-line fasta.
AllFastaStats\t| afs\tGenerate a summary of multi-line fasta statistics.
MultipleFastaStats\t| mfs\tGenerate a summary of multi-line fasta statistics (Multiple).
ConcatenateFasta\t| ccf\tMake a concatenated fasta file (Multiple).
TranslateSequence\t| tls\tFind the translated sequences as a protein and open reading frames (ORFs).

Use <module> --help for module usage.'''

module_map = {
    'Multi2Single': 'multi2single.py',
    'm2s': 'multi2single.py',
    'RenameId': 'renameid.py',
    'rid': 'renameid.py',
    'PrefixRename': 'prfxrename.py',
    'prn': 'prfxrename.py',
    'PrefixSelectRename': 'prfxselrename.py',
    'psr': 'prfxselrename.py',
    'IdExtract': 'idextract.py',
    'idx': 'idextract.py',
    'IdExtractLocation': 'idextloct.py',
    'iel': 'idextloct.py',
    'IdExtractLocationMultiple': 'idextloctmlt.py',
    'iem': 'idextloctmlt.py',
    'ReverseComplement': 'revcomplt.py',
    'rcp': 'revcomplt.py',
    'FindCountDuplication': 'findcntdupl.py',
    'fcd': 'findcntdupl.py',
    'RemoveDuplication': 'removedupl.py',
    'rvd': 'removedupl.py',
    'SubsetFasta': 'subsetfa.py',
    'ssf': 'subsetfa.py',
    'ExtractPattern': 'extractptrn.py',
    'xpt': 'extractptrn.py',
    'EachFastaStats': 'eachfastats.py',
    'efs': 'eachfastats.py',
    'AllFastaStats': 'allfastats.py',
    'afs': 'allfastats.py',
    'MultipleFastaStats': 'asmstatsunlm.py',
    'mfs': 'asmstatsunlm.py',
    'ConcatenateFasta': 'concatenate.py',
    'ccf': 'concatenate.py',
    'TranslateSequence': 'translatedna.py',
    'tls': 'translatedna.py'
}

if __name__ == '__main__':
    if len(sys.argv) < 2 or sys.argv[1] in ['--h', '--help']:
        print(usage)
        sys.exit(0)

    module = sys.argv[1]
    if module not in module_map:
        print('Unexpected module. Use --h for help.')
        sys.exit(0)

    script_path = f'{sys.path[0]}/scripts/{module_map[module]}'
    parameters = ' '.join(sys.argv[2:])
    subprocess.run(f'python3 {script_path} {parameters}', shell=True)

