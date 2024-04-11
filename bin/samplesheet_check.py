#!/usr/bin/env python3

import sys
import shutil

samplesheet = sys.argv[1]
validated_samplesheet = samplesheet.replace('.csv', '.valid.csv')
shutil.copy(samplesheet, "samplesheet.valid.csv")