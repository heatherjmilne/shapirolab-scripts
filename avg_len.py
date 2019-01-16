#!/usr/bin/python

import sys, os
# Used to calculate the average length of fragments from the lgdistribution.txt file produced by mapdamage.

def calc_avg_length(filename):
	len_sum = 0
	count = 0
	vals = []
	histogram = open(filename, "r")
	for line in histogram:
		if line[0] == "+" or line[0] == "-":
			strand, length, freq = line.strip("\n").split("\t")
			len_sum += float(length)*float(freq)
			count += float(freq)
			vals.append((float(length), float(freq)))
	avg = round(len_sum/count, 3)
	print "Average: ", avg
	print "Count: ", count
	return avg, vals

def calc_std_dev(avg, vals):
	#vals should be a (value, frequency) tuple.
	stdev = 0
	count = 0
	if len(vals) < 2:
		stdev = 0
	else:
		for pair in vals:
			value, freq = pair[0], pair[1]
			count += freq
			stdev = freq*((value - avg)**2)
		stdev = (stdev/(count - 1))**(0.5)
	print "Std Dev: ", stdev
	return stdev

def write_output(avg, stdev, filename):
	output = open(filename, "w+")
	output.write("Average length: \t")
	output.write(str(avg))
	output.write("\n")
	output.write("Std Dev: \t")
	output.write(str(stdev))
	output.write("\n")
	output.close()
	print "Output written to ", filename

i=1
while i < len(sys.argv):
	inputfile = sys.argv[i]
	pathname = os.path.abspath(os.path.dirname(inputfile))
	#print(pathname)
	outputfile = pathname + "/avg_length.txt"
	#print(outputfile)
	hist_avg, hist_values = calc_avg_length(inputfile)
	hist_stdev = calc_std_dev(hist_avg, hist_values)
	write_output(hist_avg, hist_stdev, outputfile)
	i+=1




