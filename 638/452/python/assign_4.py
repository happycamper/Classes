#!/usr/bin/python

matrix_V = [1,2] 
matrix_Z = [
	[2,1],
	[1,3]
	]
matrix_I = []

det_Z = 1.0/(matrix_Z[0][0]*matrix_Z[1][1] - matrix_Z[0][1]*matrix_Z[1][0])

temp = matrix_Z[0][0]
matrix_Z[0][0] = matrix_Z[1][1]
matrix_Z[1][1] = temp

temp = matrix_Z[0][1] * -1
matrix_Z[0][1] = matrix_Z[1][0] * -1
matrix_Z[1][0] = temp

rowCount = 0
colCount = 0

#invert a 2x2 matrix
for row in matrix_Z:
	for col in row:
		col *= det_Z
		matrix_Z[rowCount][colCount] = col
		colCount += 1	
	rowCount += 1
	colCount = 0

rowCount = 0
#calculate I
for row in matrix_Z:
	temp = matrix_Z[rowCount][0]*matrix_V[0] + matrix_Z[rowCount][1]*matrix_V[1]	
	matrix_I.append(temp)
	print "I%d %f" % (rowCount,temp) 
	rowCount += 1
