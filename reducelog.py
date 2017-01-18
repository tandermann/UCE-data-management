#written by Tobias Hofmann (tobias.hofmann@bioenv.gu.se)
import os
import argparse
import csv

# Complete path function
class CompletePath(argparse.Action):
	"""give the full path of an input file/folder"""
	def __call__(self, parser, namespace, values, option_string=None):
		setattr(namespace, self.dest, os.path.abspath(os.path.expanduser(values)))

# Get arguments
def get_args():
	parser = argparse.ArgumentParser(
		description="reduce mcmc log file to the user-defined columns",
		formatter_class=argparse.ArgumentDefaultsHelpFormatter
	)
	parser.add_argument(
		'--mcmclog',
		required=True,
		action=CompletePath,
		help='The MCMC log file'
	)
	parser.add_argument(
		'--first_x',
		type=int,
		default=5,
		help='Give the number of columns you want to reduce the logfile to. Additionally you can add specific columns of interest using the --column_names flag where you can specify the names of the additional columns you want to extract'
	)
	parser.add_argument(
		'--column_names',
		type=str,
		default=False,
		help='Specify the names of the columns you want to extract additionally to the columns from the --first_x_lines command. Input the names separated with "," e.g. name1,name2,name3'
	)
	return parser.parse_args()

# Get arguments
args = get_args()

input_file = args.mcmclog
int_cols = args.first_x
columns = args.column_names
list_columns = columns.split(',')


def exclude_commented_lines(file_content):
	content = []
	for line in file_content:
		li = ''.join(line)
		if not li.startswith("#"):
			content.append(line)
	return content

def reducelog(file,number_of_columns,additional_column_list):
	infile_path = file.split("/")
	workdir = "/".join(infile_path[:-1])
	infile_stem = infile_path[-1].split(".")[0]
	outfile = "%s/%s_reducedLog.log" %(workdir,infile_stem)
	print outfile
	output = open(outfile, "wb")
	outlog=csv.writer(output, delimiter='\t')
	print "Reading mcmc log file. This can take some moments (depending on file-size)..."
	with open(file, 'r') as f:
		reader = csv.reader(f, delimiter='\t')
		reader = list(reader)
		reader = exclude_commented_lines(reader)
		header = reader[0]
		target_columns = header[:number_of_columns] + additional_column_list
		body = reader[1:]
		id_list = []
		for column in header:
			if column in target_columns:
				id_list.append(header.index(column))
		out_header = []
		for element in id_list:
			out_header.append(header[element])
		outlog.writerow(out_header)
		for line in body:
			out_line = []
			for element in id_list:
				out_line.append(line[element])
			outlog.writerow(out_line)





reducelog(input_file,int_cols,list_columns)