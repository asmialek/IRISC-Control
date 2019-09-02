import os
import sys
from decimal import *

getcontext().prec = 9

filename = sys.argv[1]

with open(filename, 'r') as f:
	file_lines = f.readlines()

new_file = ['Time,X,Y,Z\n']

first_line = file_lines[0].split(',')
first_time = Decimal(first_line[0])

i = float(0.00)
for line in file_lines:
	print(line)
	temp_line = line.split(',')
	new_time = Decimal(temp_line[0]) - first_time
	new_line = str(i) + ',' + ','.join(temp_line[1:])
#	new_line = temp_line[1] + '\n'
	new_file.append(new_line)
	print(new_line)
	i += 0.01

new_filename = filename.split('.')[0] + '_new.csv'
with open(new_filename, 'w') as f:
	f.writelines(new_file)

