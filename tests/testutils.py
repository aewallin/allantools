# read a simple data-file with phase or frequency numbers on each line
def read_datafile(filename):
	p=[]
	with open(filename) as f:
		for line in f:
			if line.startswith("#"): # skip comments
				pass
			else:
				p.append( float(line) )
	return p

# read a result-file, produced by copy/paste from Stable32
# note: header-lines need to be manually commented-out with "#"
def read_resultfile(filename):
	rows=[]
	row=[]
	with open(filename) as f:
		for line in f:
			if line.startswith("#"):
				pass
			else:
				row = []
				l2 = line.split(" ")
				l2 = filter(None, l2)
				for n in range(len(l2)):
					row.append( float(l2[n]) )
				rows.append(row)
	return rows

# test a deviation function by:
# - running the function on the datafile
# - reading the correct answers from the resultfile
# - checking that tau, n, and dev are correct 
def test( function, datafile, datainterval, resultfile, verbose=0, tolerance = 1e-4):
	# if Stable32 results were given with more digits we could decrease tolerance
	
	phase = read_datafile(datafile)
	print "Read ",len(phase)," phase values from ", datafile
	devresults = read_resultfile(resultfile)
	print "Read ",len(devresults)," rows from ", resultfile
	taus = []
	devs = []
	ns   = []

	# parse textfile produced by Stable32    
	for row in devresults:
		if len(row)==7: # typical ADEV result file has 7 columns of data
			tau_n = row[0] # tau in number of datapoints
			tau_s = row[1] # tau in seconds
			taus.append( tau_s )
			n = row[2] # n averages
			a = row[5] # deviation
			# Note we don't read the error-columns, since they are mostly zero for 'all tau' runs of Stable32   
			devs.append(a)
			ns.append(n)
		elif len(row)==4: # the MTIE/TIErms results are formatted slightly differently
			tau_n = row[0] # tau in number of datapoints
			tau_s = row[1] # tau in seconds
			taus.append( tau_s )
			n = row[2] # n averages
			a = row[3] # MTIE or TIErms
			devs.append(a)
			ns.append(n)
	# run allantools algorithm   
	(taus2,devs2,errs2,ns2) = function(phase, datainterval, taus)
	# check that allantools and Stable32 agree on length of DEV, Tau, and N results
	assert( len(taus) == len(taus2) )
	assert( len(devs) == len(devs2) )
	assert( len(ns) == len(ns2) )
	n_errors = 0 # number of errors/problems detected
	for (t1,a1,n1,t2,a2,n2) in zip(taus,devs,ns, taus2,devs2,ns2):
		# Check that allantools and Stable32 give exactly the same Tau and N results
		try:
			assert( t1 == t2 )
		except:
			print "ERROR t1=", t1, " t2=", t2
			n_errors = n_errors + 1
		try:
			assert( n1 == n2 )
		except:
			n_errors = n_errors + 1
			print "ERROR n1=", n1, " n2=",n2," at t1=", t1, " t2=", t2
		
		# check the DEV result, with a given relative error tolerance        
		rel_error = (a2-a1)/a1
		#tolerance = 1e-4 # if Stable32 results were given with more digits we could decrease tol
		try:
			assert( abs(rel_error) < tolerance )
			if verbose:
				print "OK %d %d  \t %0.6f \t %0.6f \t %0.6f" % (t1,n1,a1,a2, rel_error)
		except:
			n_errors = n_errors + 1
			print "ERROR %d  %d %0.6f \t %0.6f \t %0.6f" % (t1,n1,a1,a2, rel_error)
			n_errors = n_errors + 1
	assert( n_errors == 0) # no errors allowed!
	print "test of function ",function, " Done. Errors=", n_errors
	
if __name__ == "__main__":
	pass
