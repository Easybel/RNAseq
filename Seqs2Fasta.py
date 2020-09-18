import sys

# this script:
#    - takes a list of sequences
#    - and prints it by adding lines inbetween that fulfill the fasta convention

#  input argument: name of the file
print(sys.argv)

filename_in = sys.argv[1] + ".txt"
filename_out = sys.argv[1] + ".fasta"

with open(filename_in, 'r+') as f_in:
	with open(filename_out, 'w+') as f_out:
		i = 1
		for line in f_in.readlines():
			#f_out.write(f'>{i:2d}\n')
			f_out.write('>' + str(i) + ' overrepresented_sequence' + '\n')
			f_out.write(line)
			i += 1
