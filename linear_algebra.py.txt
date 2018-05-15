from __future__ import absolute_import 
import csv
import linear_algebra

class linear_algebra_error(Exception):
    '''
    This class is used for defining exception for linear_algebra
    '''
    pass

class linear_algebra:
    '''
    This class is used for defining library functions for Linear Algebra
    Including Matrix Operation, Basic statistical operations, Solving Linear Equations
    and Principal Component Analysis
    '''
    
    ENABLE_LOG = False
    @staticmethod
    def log(message):
        '''
        This method is used for logging purpose. Ny defualt log is turned off, to turn on
        Change ENABLE_LOG to True
        '''
        if linear_algebra.ENABLE_LOG:
            print message    
        
        
   
    @staticmethod    
    def get_matrix_dimensions(matrix):
        '''  
        This method returns dimension of a matrix as a  tuple (rows,cols) as output
        '''
        return len(matrix), len(matrix[0]) 
        
    @staticmethod
    def print_matrix(matrix):
        '''  
        This method prints matrix in a fancy format
        '''
        n = len(matrix)
        m = len(matrix[0])
        for i in range(0, n):
            line = ""
            for j in range(0, m):
                line += str(matrix[i][j]) + "\t"
                if j == n-1:
                    line += "| "
            print(line)
        print("")
    
    @staticmethod
    def get_matrix_transpose(matrix):
        '''  
        This method returns transpose of a matrix
        '''
        result = [[matrix[j][i] for j in xrange(len(matrix))] for i in xrange(len(matrix[0]))]
        return result         

    
    @staticmethod
    def get_matrix_determinant(matrix):
        '''  
        This method returns returns determinant of a input matrix
        '''
        if len(matrix) == 2:
            return matrix[0][0]*matrix[1][1] - matrix[0][1] * matrix[1][0]

        det = 0.0
        for c in range(len(matrix)):
            det += ((-1.0)**c)*matrix[0][c]* linear_algebra.get_matrix_determinant(linear_algebra.get_matrix_minor(matrix,0,c))
        return det

    @staticmethod
    def get_matrix_minor(matrix,i,j):
        '''  
        This method returns returns minors of a matrix
        '''
        return [row[:j] + row[j+1:] for row in (matrix[:i] + matrix[i+1:])] 

    @staticmethod
    def get_matrix_inverse(matrix):
        '''  
        This method returns returns inverse of of a matrix. If matrix is not invertible i.e. 
        determinant is 0, it will raise a linear_algebra_error
        '''
        det = linear_algebra.get_matrix_determinant(matrix)
        
        if (det == 0.0):
           raise linear_algebra_error("Matrix is singular, can not be inversed")
        linear_algebra.log("determinant=" + str(det))
        #special case for 2x2 matrix:
        if len(matrix) == 2:
            return [[matrix[1][1]/det, -1 * matrix[0][1]/det],
                   [-1*matrix[1][0]/det, matrix[0][0]/det]]

    #find matrix of cofactors
        cofactors = []
        for r in range(len(matrix)):
            cofactorRow = []
            for c in range(len(matrix)):
                minor = linear_algebra.get_matrix_minor(matrix,r,c)
                cofactorRow.append(((-1.0)**(r+c)) * linear_algebra.get_matrix_determinant(minor))
            cofactors.append(cofactorRow)
        cofactors = linear_algebra.get_matrix_transpose(cofactors)
    
        for r in range(len(cofactors)):
            for c in range(len(cofactors)):       
                cofactors[r][c] = cofactors[r][c]/det
        return cofactors
        
    @staticmethod
    def get_matrix_norm(matrix, factor):
        '''
        This method returns norm of a matrix
        '''
        n = sum([ matrix[i][j] ** factor for i in range(len(matrix)) for j in range(len(matrix[0]))])
        return n ** (1.0/ factor) 

    @staticmethod
    def get_matrix_sum(matrix1, matrix2):
        '''
        This method does summation of two matrices matrix1, matrix2
        '''
        row_count1 = len(matrix1)
        row_count2 = len(matrix2)
        column_count1 = len(matrix1[0])
        column_count2 = len(matrix2[0])
        if (row_count1 != row_count2):
            raise linear_algebra_error("Row counts of matrices do not match: Sum operation can not be performed")
            
        if (column_count1 != column_count2):
            raise linear_algebra_error("Column counts of matrices do not match: Sum operation can not be performed")
        return [[matrix1[i][j] + matrix2[i][j]  for j in range(len(matrix1[0]))] for i in range(len(matrix1))] 
    
    @staticmethod
    def get_matrix_diff(matrix1, matrix2):
        '''
        This method does summation of two matrices matrix1, matrix2
        '''
        row_count1 = len(matrix1)
        row_count2 = len(matrix2)
        column_count1 = len(matrix1[0])
        column_count2 = len(matrix2[0])
        if (row_count1 != row_count2):
            raise linear_algebra_error("Row counts of matrices do not match: Sum operation can not be performed")
            
        if (column_count1 != column_count2):
            raise linear_algebra_error("Column counts of matrices do not match: Sum operation can not be performed")
            
        return [[matrix1[i][j] - matrix2[i][j]  for j in range(len(matrix1[0]))] for i in range(len(matrix1))] 
        
    @staticmethod
    def get_matrix_product(matrix1, matrix2):
        '''
        This method returns product of two matrices matrix1, matrix2
        '''
        
        row_count1 = len(matrix1)
        row_count2 = len(matrix2)
        column_count1 = len(matrix1[0])
        column_count2 = len(matrix2[0])  
        if (column_count1 != row_count2 ):
            raise linear_algebra_error("Column count of first matrix does not match row count of second matrix: Matrix product operation can not be performed") 
            
        result = [[sum(m1*m2 for m1,m2 in zip(matrix1_row, matrix2_column)) for matrix2_column in zip(*matrix2)] for matrix1_row in matrix1]
        return result
        
    @staticmethod
    def get_matrix_identity(length):
        '''
        This method returns an indentity matrix size length * length
        '''
        return  [[1 if i==j else 0 for j in range(length)] for i in range(length)] 
        
    @staticmethod
    def get_matrix_scalar_product(val, matrix):
        '''
        This method does scalar product val of elements of matrix
        '''
        return  [[val * matrix[i][j] for j in range(len(matrix[0]))] for i in range(len(matrix))]  
    
    @staticmethod
    def removeZero(element):
        return element[0] <> 0.0  

    @staticmethod
    def get_eigen_values(matrix):
        '''
        This method returns a list of eigen values of a given matrix
        Eigen value is calculate using QR algorithm
        '''
        
        A = matrix
        R = A
        QT = [[]]
    
        iteration = 1
        tol =0.01
        while abs(tol) >= 0.00001 and iteration <= 100:
            linear_algebra.log(" Working for iteration=" + str(iteration))
            R = A
            QT = [[]]
            for k in xrange(len(R[0]) - 1):
                column = linear_algebra.get_matrix_transpose([list(zip(*R)[k])])
                linear_algebra.log("column[" + str(k) + "] =" + str(column))
                nm = linear_algebra.get_matrix_norm(column, 2)
                linear_algebra.log("nm=" + str(nm))
                if nm == 0.0:
                    std_column = column
                else:
                    std_column = linear_algebra.get_matrix_scalar_product(1.0/nm, column)       
                linear_algebra.log( "std_column[" + str(k) + "] =" + str(std_column))
                d = linear_algebra.get_matrix_norm(std_column[k:][:], 2)
                if (std_column[k][0] > 0.0) :
                    d = -1.0 * d
                linear_algebra.log( "d=" + str(d))
                vk = ( 0.5 * (1-std_column[k][0]/d)) ** 0.5
                t = -d * vk
                linear_algebra.log("t=" + str(t))
                v = [[0.0 if i <k else vk if i==k else std_column[i][0] / (2 * t)  ] for i in xrange(len(std_column))]
                linear_algebra.log("v=" + str(v))
                identity_matrix = linear_algebra.get_matrix_identity(len(std_column))
                vtrans = linear_algebra.get_matrix_transpose(v)
                linear_algebra.log("vtrans=" + str(vtrans))
                multi_v_vtrans = linear_algebra.get_matrix_product(v, vtrans)
                p = linear_algebra.get_matrix_diff(identity_matrix, linear_algebra.get_matrix_scalar_product(2.0, multi_v_vtrans))
                linear_algebra.log("vtrans=" + str(vtrans))
                R = linear_algebra.get_matrix_product(p, R)
                linear_algebra.log("R=" + str(R))
                if k==0:
                    QT = p
                else:
                    QT = linear_algebra.get_matrix_product(p, QT)
                linear_algebra.log( "QT=" + str(QT))
            
            Q = linear_algebra.get_matrix_transpose(QT)
            linear_algebra.log("Q=" + str(Q))
            linear_algebra.log("R=" + str(R))
            res = linear_algebra.get_matrix_product(Q, R)
            linear_algebra.log("res=" + str(res))
            newA = linear_algebra.get_matrix_product(R, Q)
            linear_algebra.log("newA=" + str(newA))
            diff_matrix = linear_algebra.get_matrix_diff (newA, A)
            tol = max([abs(diff_matrix[i][j]) if i==j else 0 for i in xrange(len(diff_matrix)) for j in xrange(len(diff_matrix[0]))])
            A = newA
            linear_algebra.log( " After iteration=" + str(iteration) + " diff_matrix = " + str(diff_matrix) + " tolerance=" + str(tol))
            iteration += 1
    
        linear_algebra.log( "Final A=" + str(A))
        eigen_value_list = [A[k][k] for k in xrange(len(A))] 
        linear_algebra.log("eigen_value_list=" + str(eigen_value_list)) 
        return eigen_value_list 
    
    @staticmethod
    def get_eigen_vector(matrix, eigen_value):
        '''
        This method finds eigen vector for a given eigen value eigen_value
        This is calculated using Inverse Iteration https://en.wikipedia.org/wiki/Inverse_iteration
        '''
        
        eigen_value += 0.000001
        length = len(matrix)
        #bprev_vector = [[1 for i in xrange(length)]]
        bprev_vector = [[1] for i in xrange(length)]
    
        linear_algebra.log( "bprev_vector=" + str(bprev_vector))
   
        identity_matrix = linear_algebra.get_matrix_identity(length)
        m = linear_algebra.get_matrix_scalar_product(eigen_value, identity_matrix)
        linear_algebra.log( "matrix=" + str(matrix))
        linear_algebra.log( "m=" + str(m))
    
        s = linear_algebra.get_matrix_diff(matrix, m)
        linear_algebra.log( "s=" + str(s))
        inv = linear_algebra.get_matrix_inverse(s)
        nrm = linear_algebra.get_matrix_norm(inv, 2)
        inv_normalized = inv
        linear_algebra.log( "inv_normalized=" + str(inv_normalized))
        if (nrm > 0.0):
            inv_normalized = linear_algebra.get_matrix_scalar_product(1.0/nrm, inv)
         
        iteration = 1
        tol =0.01
        while abs(tol) >= 0.00001 and iteration <= 100:
            b_vector = linear_algebra.get_matrix_product(inv_normalized, bprev_vector)
            linear_algebra.log( " After iteration=" + str(iteration) + " b_vector=" + str(b_vector))
            diff_vector = linear_algebra.get_matrix_diff (b_vector, bprev_vector)
            tol = max([abs(diff_vector[i][0]) for i in xrange(len(diff_vector))])
            linear_algebra.log( " After iteration=" + str(iteration) + ", diff_vector=" + str(diff_vector) + " tolerance=" + str(tol))
            iteration += 1
            bprev_vector = b_vector
        
        abs_vector =   [bprev_vector[i][0] for i in xrange(len(bprev_vector))]
        linear_algebra.log( "abs_vector=" + str(abs_vector))
        max_val = max(abs_vector)  
        linear_algebra.log( "max_val=" + str(max_val))    
        bprev_vector_zero = [[bprev_vector[i][0] if max_val/abs(bprev_vector[i][0]) <100000 else 0.0] for i in xrange(len(bprev_vector))]    
        linear_algebra.log( "bprev_vector_zero=" + str(bprev_vector_zero) ) 
        bprev_vector_no_zero = filter(linear_algebra.removeZero,  bprev_vector_zero)
        linear_algebra.log( "bprev_vector_no_zero=" + str(bprev_vector_no_zero) )
        min_val = 1.0
        if (len(bprev_vector_no_zero) > 0) :
            min_val = min([abs(bprev_vector_no_zero[i][0]) for i in xrange(len(bprev_vector_no_zero))])
        linear_algebra.log( "min_val=" + str(min_val) )
        final_vector =  [[bprev_vector_zero[i][0]/min_val] for i in xrange(len(bprev_vector_zero))]
        less_than_zero_vector = [final_vector[i][0] < 0 for i in xrange(len(final_vector))]
        linear_algebra.log( "less_than_zero_vector=" + str(less_than_zero_vector))
        less_than_zero_vector_length = sum(less_than_zero_vector)
        linear_algebra.log( "less_than_zero_vector=" + str(less_than_zero_vector_length))
        if less_than_zero_vector_length == len(final_vector):
            final_vector =  [[-1.0 * final_vector[i][0]] for i in xrange(len(bprev_vector_zero))]
        linear_algebra.log( "final_vector=" + str(final_vector))
    
        return final_vector 
        
    @staticmethod
    def swap_rows(matrix, i, j):
        '''
        This method is used for swapping two rows i, j of a matrix
        '''
        matrix[i], matrix[j] = matrix[j], matrix[i]
        return matrix

    @staticmethod
    def get_matrix_rank(matrix):
        ''' 
        This method finds rank of a matrix
        '''
        no_of_columns = len(matrix[0])
        no_of_rows = len(matrix)
        rank = min(no_of_rows, no_of_columns)
    
        if no_of_rows > no_of_columns:
            matrix = linear_algebra.get_matrix_transpose(matrix)
            no_of_columns, no_of_rows = no_of_rows, no_of_columns

        for row in xrange(rank):
            if matrix[row][row] != 0.0:
                for column in xrange(no_of_columns):
                    if (column != row):               
                        multiplication_factor = matrix[column][row] /matrix[row][row];
                        for i in xrange(rank):
                            matrix[column][i] -= multiplication_factor * matrix[row][i];
            else:
                reduce = True
                for i in range(row + 1,no_of_rows):
                    if matrix[i][row] == 0.0:
                        linear_algebra.swap_rows(matrix, row, i);
                        reduce = False;				
                        break

                if reduce:
                    for i in xrange(no_of_rows):
                        matrix[i][row],matrix[i][rank -1] = matrix[i][rank -1],matrix[i][row]
                no_of_rows -= 1


	    count_zero_elements_rows = 0
	    for i in matrix:
		    if i==[0]* no_of_columns:
			    count_zero_elements_rows += 1

        return (rank - count_zero_elements_rows)
        
    @staticmethod    
    def solve_gauss_jordan(A):
        '''
        This function solves a set of linear equations using Gauss Jordan method 
        in which a given augmented matrix A is converted to row echelon form
        '''
        
        n = len(A)
        #rank = linear_algebra.get_matrix_rank(A)
        #if (rank < n):
        #    raise linear_algebra_error("Rank of augmented matrix: " + str(rank) + " less than number of rows. Unique solution is not possible" )

        for i in range(0, n):
            # Search for maximum in this column
            maxEl = abs(A[i][i])
            maxRow = i
            for k in range(i+1, n):
                if abs(A[k][i]) > maxEl:
                    maxEl = abs(A[k][i])
                    maxRow = k

            # Swap maximum row with current row (column by column)
            for k in range(i, n+1):
                tmp = A[maxRow][k]
                A[maxRow][k] = A[i][k]
                A[i][k] = tmp

            # Make all rows below this one 0 in current column
            for k in range(i+1, n):
                c = -A[k][i]/A[i][i]
                for j in range(i, n+1):
                    if i == j:
                        A[k][j] = 0
                    else:
                        A[k][j] += c * A[i][j]

        # Solve equation Ax=b for an upper triangular matrix A
        x = [0 for i in range(n)]
        for i in range(n-1, -1, -1):
            x[i] = A[i][n]/A[i][i]
            for k in range(i-1, -1, -1):
                A[k][n] -= A[k][i] * x[i]
        return x
        
    
    @staticmethod
    def gsiedel(A,parm):
        n=len(A)
        L = [[0 for j in range(n)] for i in range(n)]   
        U = [[0 for j in range(n)] for i in range(n)]   
        B = [0 for i in range (n)]   
    
        for i in range(n):
            for j in range(n + 1):            
                if (j == n):
                    B[i] = A[i][j]
                else:
                    if (j <= i): 
                        L[i][j] = A[i][j]
                    else:                    
                        U[i][j] = A[i][j]
   
        if parm == 1:
            return L
        else:
            if parm == 2:
                return U
            else:
                return B 
                
    @staticmethod            
    def compare(a,b):
        compareresult=0
   
        for i in range(len(a)):   
            ai = a[i][0]
            ai= round(ai,2)
            if (round(a[i][0],3) == round(b[i][0],3)):
                compareresult = 1
            else:
                compareresult = 0
                break
        if compareresult == 1:
            return 1
        else:
            return 0  
            
    @staticmethod            
    def solveeq(T,C,X):
        ST = X
        brk=0
        count=0
        linear_algebra.log( "count=" + str(count) + " brk =" + str(brk) + " count =" + str(count))
    
        while brk == 0:        
            ix = linear_algebra.get_matrix_product(T, X)
            X =  linear_algebra.get_matrix_sum(ix,C)
            compres= linear_algebra.compare (X,ST) 
            linear_algebra.log( "count=" + str(count) + " X=" + str(X) + ": ST=" + str(ST) + ": compres=" + str(compres))
            if compres == 1:
                brk=1
                break
            else:
                ST=X
                count += 1          
    
        return X
    
    @staticmethod
    def solve_gauss_seidel(A): 
        '''
        This function solves a set of linear equations using Gauss Jordan method 
        in which a given augmented matrix A is converted to row echelon form
        '''   
        
        n = len(A)
        #rank = linear_algebra.get_matrix_rank(A)
        #if (rank < n):
        #    raise linear_algebra_error("Rank of augmented matrix: " + str(rank) + " less than number of rows. Unique solution is not possible" )
        
        # This code gets the Lower Matrix of the Augmented Matrix##
        linear_algebra.log( "A=" + str(A))
        parm= 1
        L= linear_algebra.gsiedel(A,parm)
        linear_algebra.log( "L=" + str(L))

        # This code gets the Strictly upper Matrix of the Augmented Matrix#
        parm=2
        U = linear_algebra.gsiedel(A,parm)
        linear_algebra.log( "U=" + str(U))
    
        # This code gets the Value Matrix of the Augmented Matrix#    
        parm=3
        B = linear_algebra.gsiedel(A,parm)
        B= list(map(list, zip(B)))
        linear_algebra.log( "B=" + str(B))
    
        # The Following Lines of Code will Call Matrix inverse Function and get the Inverse of Matrix L#    
        Linv = linear_algebra.get_matrix_inverse(L)  
        linear_algebra.log( "Linv=" + str(Linv))        
        # The Following Code will Call Matrix Multiplication Code and Get the Resultant Matrix #
    
        HT = linear_algebra.get_matrix_product(Linv, U)
        s = -1
        T = linear_algebra.get_matrix_scalar_product(s,HT)
        linear_algebra.log( "T=" + str(T))
        # The following lines of code will calculate the C matrix #
    
        C = linear_algebra.get_matrix_product(Linv, B)
        linear_algebra.log( "C=" + str(C))
      
        # The Following Lines of Code will Create the X Matrix ###############
    
        HX = [1.0 for i in xrange(len(A))]
        X = list(map(list, zip(HX))) 
        linear_algebra.log( "X=" + str(X))
        
        # The Following Lines of code will :
        #1) Call a function which will take the T Matrix, C Matrix and Initial X matrix
        #2) Calculate the TX+C matrix .
        #3) Compare the result with the Previous Value 
        #4) find out the Value of unknowns 

        Result = linear_algebra.solveeq(T,C,X)
        return Result
            
    @staticmethod    
    def pca(matrix, dimensions):
        '''
        This method is used for PCA transformation of data
        Input matrix used for the data
        Input dimensions used for how may dimensions pca should be transformed
        '''
        features = []
        for k in xrange(len(matrix[0])):
            feature = list(zip(*matrix)[k])
            mn = linear_algebra.get_mean(feature)
            mean_diff_feature = [feature[i] - mn for i in xrange(len(feature))]
            features.append(mean_diff_feature)        
        linear_algebra.log( "features="  + str(features))
        covariance_matrix = []
        for i in xrange(len(features)):
            row_covariance = []
            for j in xrange(len(features)):
                if (i==j):
                    row_covariance.append(linear_algebra.get_variance(features[i]))
                else :
                    row_covariance.append( linear_algebra.get_covariance(features[i], features[j]))
            covariance_matrix.append(row_covariance)
    
        
        linear_algebra.log( "covariance_matrix="  + str(covariance_matrix))
        eigen_value_list = linear_algebra.get_eigen_values(covariance_matrix)
        linear_algebra.log( "eigen_value_list="  + str(eigen_value_list))
        eigen_value_list.sort(reverse = True)
        linear_algebra.log( "descending sorted eigen_value_list="  + str(eigen_value_list)) 
        #eigen_vector_list = [linear_algebra.get_eigen_vector(covariance_matrix, eigen_value_list[i]) for i in xrange(len( eigen_value_list))]
        eigen_vector_list_all = [linear_algebra.get_eigen_vector(covariance_matrix, eigen_value_list[i]) for i in xrange(len( eigen_value_list))]
        linear_algebra.log(" I am here")
        linear_algebra.log( "eigen_vector_list_all=" + str(eigen_vector_list_all))
        #eigen_vector_list = [linear_algebra.get_eigen_vector(covariance_matrix, eigen_value_list[i]) for i in xrange(len( eigen_value_list))][0]
        #linear_algebra.log( "descending eigen vector list="  + str(eigen_vector_list))
        #filtered_eigen_vector_list  = eigen_vector_list[0][:dimensions]
        #filtered_eigen_vector_list = eigen_vector_list[:dimensions]
        #filtered_eigen_vector_list = eigen_vector_list[:]
        filtered_eigen_vector_list = [[eigen_vector_list_all[i][j][0] for j in xrange(len(eigen_vector_list_all[0]))] for i in xrange(dimensions)]
        linear_algebra.log( "filtered_eigen_vector_list="  + str(filtered_eigen_vector_list))
        row_feature_vector = filtered_eigen_vector_list
        #row_feature_vector =  linear_algebra.get_matrix_transpose(filtered_eigen_vector_list)
        linear_algebra.log( "row_feature_vector="  + str(row_feature_vector)) 
        #row_data_adjust = get_matrix_transpose(features)
        row_data_adjust = features
        linear_algebra.log( "row_data_adjust="  + str(row_data_adjust))
        final_data = linear_algebra.get_matrix_transpose(linear_algebra.get_matrix_product(row_feature_vector, row_data_adjust)) 
        linear_algebra.log( "final_data="  + str(final_data))
        return final_data        
        
    @staticmethod    
    def convert_matrix_float_type(matrix):
        return [[float(matrix[i][j]) for j in range(len(matrix[0]))] for i in range(len(matrix))] 
    
    @staticmethod    
    def load_from_csv(file_name, header = False):
        try:
            data = list(csv.reader(open(file_name))) 
        except:
            raise linear_algebra_error("Error in reading file: " + file_name)
        
        data_filter = data
        if header == True:
            data_filter = data[1:][:]
        return data_filter
        
    @staticmethod
    def get_mean(arr):
        return sum(arr)/len(arr)

    @staticmethod
    def get_variance(arr):
        mean_value = linear_algebra.get_mean(arr)
        x_diff_mean_sqr = [(x - mean_value) ** 2 for x in arr]
        return sum(x_diff_mean_sqr)/(len(arr) -1)

    @staticmethod
    def get_covariance(arr1, arr2):
        mean1 = linear_algebra.get_mean(arr1)
        mean2 = linear_algebra.get_mean(arr2)
        product_diff_mean1_diff_mean2 = [(x-mean1) * (y-mean2) for x,y  in zip(arr1,  arr2)]
        return sum(product_diff_mean1_diff_mean2)/(len(arr1) -1)        

  

if __name__ == "__main__":

    #multiple_linear_regression()
    #matrix = [[1,2], [3,4]]
    matrix = [[4.0,7.0],[2.0,6.0]]
    inv_matrix = linear_algebra.get_matrix_inverse(matrix)
    dim = linear_algebra.get_matrix_dimensions(matrix)
    print str(dim)
    print str(inv_matrix)
    #n = norm(matrix)
    #print "norm=" + str(n)
   #matrix = [[2.0,1.0,1.0],[1.0,2.0,1.0],[1.0,1.0,2.0]]
   # For eigen value
    matrix = [[6.0,3.0,-4.0],[-5.0,-2.0,2.0],[0.0,0.0,-1.0]]
   #get_eigen_value(matrix)
   #x = get_matrix_transpose(matrix)
   #print "x=" + str(x)
    #matrix = [[1.0,-2.0,-6.0],[-2.0,2.0,-5.0],[2.0,1.0,8.0]]
    eigen_values = linear_algebra.get_eigen_values(matrix)
    print "eigen_values=" + str(eigen_values)
    for eigen in eigen_values:
        eigen_vector = linear_algebra.get_eigen_vector(matrix, eigen)
        print "for eigen_value=" + str(eigen) + " eigen_vector=" + str(eigen_vector)
    #eigen_vector = get_eigen_vector(matrix, 3.000001)
   #eigen_vector = get_eigen_vector(matrix, 1.000001)
   #eigen_vector = get_eigen_vector(matrix, -1.000001)
   #eigen_vector = get_eigen_vector(matrix, 3.000001)
   #matrix = [[1,3],[2,2]]
   #eigen_vector = get_eigen_vector(matrix, 4.000001)
   #print "final eigen vector =" + str(eigen_vector)
   #print "Before calling product"
   #matrix1 = [[1,-1,2], [3,1,4]]
   #matrix2 = [[2,-1,3],[5,1,2],[4,6,-2]]
   #m = get_matrix_product(matrix1, matrix2)
   #m = get_matrix_product(matrix2, matrix2)
   #print "Multiply=" + str(m)
   # For loading data
   #load_csv("customers.csv")
    matrix = [[2.5,2.4], [0.5, 0.7], [2.2,2.9], [1.9,2.2],[3.1,3.0], [2.3,2.7], [2.0,1.6], [1.0,1.1], [1.5, 1.6], [1.1, 0.9] ]
    linear_algebra.pca(matrix, 2)
    
    matrix_rank = [[10,   20,   10],[20,   40,   20],[30,   50,   0]]
    print "matrix=" + str(matrix_rank)
   #swap_rows(matrix_rank, 0, 1)
   #print "after swapping matrix=" + str(matrix_rank)
  # rank = get_rank(matrix_rank)
   #print "rank=" + str(rank)
    rank = linear_algebra.get_matrix_rank(matrix_rank)
    print "rank=" + str(rank)
   
    matrix2 = [[10,   20,   10], [-20,  -30,   10], [30,   50,   0]];
    rank = linear_algebra.get_matrix_rank(matrix2)
    print "matrix2 rank=" + str(rank)
   
    matrix3 = [[10,   20,   10], [20,  40,   20], [30,   60,   30]];
    rank = linear_algebra.get_matrix_rank(matrix3)
    print "matrix3 rank=" + str(rank)
   
    #matrix1 = [[1,1], [2,2]]
    #la = linear_algebra(matrix1)
    #dim = la.dimensions()
    #print "dimensions=" + str(dim)
    #la.printMatrix()
   # A= [[10,-1,2,0,6],[-1,11,-1,3,25],[2,-1,10,-1,-11],[0,3,-1,8,15]]
   # A= [[10.0,-1.0,2.0,0.0,6.0],[-1.0,11.0,-1.0,3.0,25.0],[2.0,-1.0,10.0,-1.0,-11.0],[0.0,3.0,-1.0,8.0,15.0]]
    A = [[16.0,3.0,11.0], [7.0, -11.0, 13.0]]
 #   X = linear_algebra.solve_gauss_jordan(A)
    #linear_algebra.print_matrix(X)
 #   print "X=" + str(X)
    X2 = linear_algebra.solve_gauss_seidel(A)
    print "X2=" + str(X2)
    
    Y1 = [[10,-1,2,0],[-1,11,-1,3],[2,-1,10,-1],[0,3,-1,8]]
    #Y2 =[[1], [2], [-3], [-2]]
    Y2 =[[1.0], [2.0], [-1.0], [1.0]]
    Y3 = linear_algebra.get_matrix_product(Y1, Y2)
   
    print "Y3 = " + str(Y3)
    
    Y4=[[1.349391319689485], [1.6573747353563866], [-0.918269230769231], [0.9999999999999999]]
    Y5 = linear_algebra.get_matrix_product(Y1, Y4)
    print "Y5 = " + str(Y5)
    
   
