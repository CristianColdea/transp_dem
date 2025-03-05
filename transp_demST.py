"""
Script to compute shor term transport demand.
Three methods will be coded: Furness, Fratar and a new one.
"""

# number of travels as a matrix with produced travels on lines
# and attracted travels on columns

travs = [[40, 110, 150],
         [50, 20, 30],
         [110, 30, 10]]

# the matching friction factors, same arrangement

ffs = [[0.753, 1.597, 0.753],
       [0.987, 0.753, 0.765],
       [1.597, 0.765, 0.753]]

# neutral calibration coefficients

k_ij0 = [[1, 1, 1],
         [1, 1, 1],
         [1, 1, 1]]

# auto travels cost

tca = [[0.5, 1, 1.4],
       [1.2, 0.8, 1.2],
       [1.7, 1.5, 0.7]]

# transit travels cost

tct = [[1, 1.5, 2],
       [1.8, 1.2, 1.9],
       [1.7, 1.5, 0.7]]

# auto travels duration

tda = [[3, 12, 7],
       [13, 3, 19],
       [9, 16, 4]]

# transit travels duration

tdt = [[15, 5, 12],
       [15, 6, 26],
       [20, 21, 8]]

# the future friction factors

ffs_f = [[0.753, 0.987, 1.597],
         [0.987, 0.753, 0.765],
         [1.597, 0.765, 0.753]]

# the future produced travels

P_is = [750, 580, 480]

# the future attracted travels

A_js = [722, 786, 302]

class GravitMod:
    def __init__(self, travs, ffs, k_ijs, P_is, A_js):
        self.travs = travs
        self.ffs = ffs
        self.k_ijs = k_ijs
        self.P_is = P_is
        self.A_js = A_js

    def gravmod_init(travs,ffs, k_ijs):
        """
        Method to compute gravitational model values in order to determine the
        calibration factors.
        Takes as input the travels, friction factors and neutral calibration
        coefficients matrices.
        Returns a matrix with the computed travels.
        """
    
        # check if the matrices have the same shape
        if(len(travs) != len(ffs) or (len(travs) != len(k_ijs))):
            print("The matrices doesn't match. Please fix it.")
            exit()
    
        # transpose de matrices
        travs_tt = list(zip(*travs))
        ffs_tt = list(zip(*ffs))
        travs_t = [list(sublist) for sublist in travs_tt]
        ffs_t = [list(sublist) for sublist in ffs_tt]
        # print("travs_tt, ", travs_tt)
        # print("travs_t, ", travs_t)
        # print(ffs_t)
    
        # get attracted travels sums (cycling on transposes)
        s_Aj = []   # store the attracted sums

        for item in travs_tt:
            s_Aj.append(sum(item))

        # get produced travels sums
        s_Pi = []
        for item in travs:
            s_Pi.append(sum(item))
    
        # compute travels with gravitational model
        gvals_init = []    # to store computed values
        for i in range(len(travs)):
            pdsum = 0
            for j1, j2 in zip(s_Aj, ffs[i]):
                pdsum = pdsum + j1 * j2
                # print(pdsum)
                # print(ffs[i])
            for k1 in range(len(ffs[i])):
                gvals_init.append((s_Pi[i] * ffs[i][k1] * s_Aj[k1] * k_ijs[i][k1] /
                                   pdsum))

        #print("Initial travels obtained with gravitational model, ", gvals_init)

        # check raw produced travels
        gvals_init_m0 = [gvals_init[i:i + 3] for i in range(0, len(gvals_init), 3)]
        #print("Initial travels matrix is ,", gvals_init_m0)

    
        # for p1, p2 in zip(travs, gvals_init_m0):
            #print(round(sum(p1)) == round(sum(p2)))
            #print(sum(p2))

        # round the number of travels
        gvals_init_r = []
        for item in gvals_init:
            gvals_init_r.append(round(item))
    
        #print("Rounded number of initial travels, ", gvals_init_r)

        # group flatten list 'gvals_init_r' as a matrix
        gvals_init_m = [gvals_init_r[i:i + 3] for i in range(0,
                         len(gvals_init_r), 3)]
        # print(gvals_init_m)
        print("Matrix of rounded numbers, ", gvals_init_m)

        # check produced travels sum
        # for p1, p2 in zip(travs, gvals_init_m):
            # print(sum(p1) == sum(p2))

        # check attracted travels sum
        # transpose the matrix first
        gvals_init_m_tt = list(zip(*gvals_init_m))
        # print(gvals_init_m_tt, travs_tt)
        # for a1, a2 in zip(gvals_init_m_tt, travs_tt):
            # print(sum(a1) == sum(a2))
            # print(sum(a1))
            # print(sum(a2))

        return gvals_init_m

    def furness(travs, P_is, A_js, tlr=0.01):

        """
        Method to iteratively compute the future travels distribution using
        Furness method.
        Takes as input the observed (historical) travels in the form of squared
        matrix,the future produced and attracted ones, respectively, in the
        form of one-line matrices (one for produced and one for attracted)
        and the precision (tolerance).
        Returns a matrix with adjusted travels.
        """

        print()
        print("Enter Furness method.")

        # check if the matrices have the same shape
        if(len(travs) != len(P_is) or len(travs) != len(A_js)):
            print("The matrices doesn't match. Please fix it.")
            exit()
        
        # function to compare the produced, respectively attracted travels
        # within a certain tolerance
        def comp(s_ih, s_ic, tlr):
            """
            Function within method to compare two values, within tolerance.
            Takes as inputs the lists of to be compared values
            and the precision/tolerance.
            Returns True of False.
            """
            
            # set a flag
            flag = True

            for ih, ic in zip(s_ih, s_ic):
                if(abs(ih - ic) / ih >= tlr): 
                    flag = False
                    break

            return flag
         
        cmp_flg = False  # comparison flag to govern the following cycle
        i = 0   # produced passes counter
        j = 0   # attracted passes counter
        
        # allocate initial computed travels matrix
        travsc = []

        for row in travs:
            travsc.append(row)

        while(cmp_flg == False):
                       
            # get produced travels sums on computed travels
            s_Pic = []
            for item in travsc:
                s_Pic.append(sum(item))

            print()
            print("P_is, ", P_is)
            print("s_Pic, ", s_Pic)
            
            cmp_flg = comp(s_Pic, P_is, tlr)
            print(cmp_flg)

            if (cmp_flg == False):
                ccs = []   # list to store produced travels coefficients
                for ps, pc in zip(P_is, s_Pic):
                    ccs.append(round(ps/pc, 3))

                # print("coefficients on produced travels, ", ccs)

                for x in range(len(travsc)):
                    travsc[x] = [ccs[x]*item for item in travsc[x]]
            
                i += 1

                print()
                print("travsc after pass  i = ", i, "is ", travsc)
            
            
            # *********
            # working on attracted travels

            # clear the coefficients vector
            ccs.clear()

            # transpose de matrices
            travsc_tt = list(zip(*travsc))

            # print()
            # print("travs_t, ", travs_t)
            # print("travsc_t ", travsc_t)
    
                    
            # get attracted travels sums on computed travels (cycling on transposes)
            s_Ajc = []   # store the attracted sums

            for item in travsc_tt:
                s_Ajc.append(sum(item))
            
            cmp_flg = comp(s_Ajc, A_js, tlr)
            
            if (cmp_flg == False):
                for ats, ac in zip(A_js, s_Ajc):
                    ccs.append(round(ats/ac, 3))

                #print()
                #print("travs, ", travs)
                #print("travsc, ", travsc)
                #print("coefficients on attracted travels, ", ccsj)

                for x in range(len(travsc_tt)):
                    travsc_tt[x] = [ccs[x]*item for item in travsc_tt[x]]
            
                j += 1

               #print()
               #print("travsc_tt after pass j = ", j, "is ", travsc_tt)
                

            travsc_0 = list(zip(*travsc_tt))
            travsc = [list(sublist) for sublist in travsc_0]

            #print()
            #print("travsc,  ", travsc)
            
            # update the attracted sums
                
            # get attracted travels sums on new computed travels (cycling on transposes)
            s_Ajc.clear()   # clear the computed attracted sums
            ccs.clear()    # clear the coefficients vector

            for item in travsc_tt:
                s_Ajc.append(sum(item))

            #print()
            #print("s_Ajc, ", s_Ajc)

            # update the produced sums
            s_Pic.clear()

            for item  in travsc:
                s_Pic.append(sum(item))

            #print()
            #print("s_Pic, ", s_Pic)
            cmp_flgA = comp(s_Ajc, A_js, tlr)
            print("Flag on attracted, ", cmp_flgA)
            
            cmp_flgP = comp(s_Pic, P_is, tlr)
            print("Flag on produced, ", cmp_flgP)

            if(cmp_flgA == True and cmp_flgP == True):
                cmp_flg = True
            else:
                cmp_flg = False

        
            travscr = []     # list to store rounded values, flatten form
            for item in travsc:
                for item in item:
                    travscr.append(round(item))

        #print()
        #print("Final rounded and flatten, ", travscr)
        travscrm = [travscr[i:i + 3] for i in range(0, len(travscr), 3)]
            
        #print()
        print("Final rounded matrix, Furness, ", travscrm)
        print("Historical travels matrix, ", travs)
        print("i is, ", i)
        print("j is, ", j)
        print("Exit Furness method.")
        
        return travscrm


    def fratar(travs, P_is, A_js, tlr=0.01):

        """
        Method to iteratively compute the future travels distribution using
        Fratar method.
        Takes as input the observed (historical) travels in the form of squared
        matrix,the future produced and attracted ones, respectively, in the
        form of one-line matrices (one for produced and one for attracted)
        and the precision (tolerance).
        Returns a matrix with adjusted travels.
        """

        print()
        print("Enter Fratar method.")

        # check if the matrices have the same shape
        if(len(travs) != len(P_is) or len(travs) != len(A_js)):
            print("The matrices doesn't match. Please fix it.")
            exit()
        
        # function to compare the produced, respectively attracted travels
        # within a certain tolerance
        def comp(s_ih, s_ic, tlr):
            """
            Function within method to compare two values, within tolerance.
            Takes as inputs the lists of to be compared values
            and the precision/tolerance.
            Returns True of False.
            """
            
            # set a flag
            flag = True

            for ih, ic in zip(s_ih, s_ic):
                if(abs(ih - ic) / ih >= tlr): 
                    flag = False
                    break

            return flag
         
        cmp_flg = False  # comparison flag to govern the following cycle
        
        p = 0   # passes counter
                
        # allocate initial computed travels matrix
        travsc = []

        for row in travs:
            travsc.append(row)

        # print("Current matrix is, ", travs)

        while(cmp_flg == False):
                       
            # get produced travels sums on current travels
            s_Pic = []
            for item in travsc:
                s_Pic.append(sum(item))

            # print()
            # print("P_is, ", P_is)
            # print("s_Pic, ", s_Pic)
            
            # working on attracted travels

            # clear the coefficients vector
           
            # transpose de matrix
            travsc_tt = list(zip(*travsc))
                 
            # get attracted travels sums (cycling on transposes)
            s_Ajc = []   # store the attracted sums

            for item in travsc_tt:
                s_Ajc.append(sum(item))

            cmp_flgA = comp(s_Ajc, A_js, tlr)
            # print("Flag on attracted, ", cmp_flgA)
            
            cmp_flgP = comp(s_Pic, P_is, tlr)
            # print("Flag on produced, ", cmp_flgP)

            if(cmp_flgA == True and cmp_flgP == True):
                cmp_flg = True
            else:
                cmp_flg = False

            if p == 1:
                cmp_flg = True

            if (cmp_flg == False):
                ccsi = []   # list to store produced travels growth factors
                for ps, pc in zip(P_is, s_Pic):
                    ccsi.append(round(ps/pc, 3))

                # print("coefficients on produced travels, ", ccsi)
                ccsj = []    # list to store attracted travels growth factors
                for ats, ac in zip(A_js, s_Ajc):
                    ccsj.append(round(ats/ac, 3))

                print("Current matrix is, ", travsc)
                
                # matrix to store t_ij travels
                travs_ij = []

                # matrix to store t_ji travels
                travs_ji = []

                for r in range(len(travsc)):
                    print("row is, ", travsc[r])
                    print("col is, ", travsc_tt[r])
                    asum = 0
                    psum = 0

                    # computing the sum of product between attracted and ccsj
                    for t in range(len(travsc[r])):
                        print("trav is, ", travsc[r][t])
                        print("ccsj[t], ", ccsj[t])
                        asum = asum + travsc[r][t] * ccsj[t]
                        print("asum is, ", asum)
                       
                    # computing the t_ij travels with Fratar formula
                    for t in range(len(travsc[r])):
                        t_ij = travsc[r][t] * P_is[r] * ccsj[t] / asum
                        travs_ij.append(t_ij)

                    # computing the sum of product between produced and ccsi
                    for tt in range(len(travsc_tt[r])):
                        print("trav is, ", travsc_tt[r][tt])
                        print("ccsi[tt], ", ccsi[tt])
                        psum = psum + travsc_tt[r][tt] * ccsi[tt]
                        print("psum is, ", psum)
                    # computing the t_ji travels with Fratar formula 
                    for tt in range(len(travsc_tt[r])):
                        t_ji = travsc[r][tt] * A_js[r] * ccsi[tt] / psum
                        travs_ji.append(t_ji)

                print("travs_ij is, ", travs_ij)
                print("travs_ji is, ", travs_ji)

                travsc0 = []

                for ij, ji in zip(travs_ij, travs_ji):
                    travsc0.append((ij+ji)/2)
                print("travsc0 is, ", travsc0)

                travsc.clear()
                travsc = [travsc0[x:x + 3] for x in range(0, len(travsc0), 3)]

                
                p += 1

            print()
            print("travsc after pass = ", p, "is ", travsc)
            
            
            
                                

            # travsc_0 = list(zip(*travsc_tt))
            # travsc = [list(sublist) for sublist in travsc_0]

        # print("Final rounded matrix, Fratar, ", travscrm)
        print("Historical travels matrix, ", travs)
        print("p is, ", p)
        print("Exit Fratar method.")
        
        return travsc

    def iter_wgt_dmd(travs, P_is, A_js, tlr=0.01):

        """
        Method to iteratively weight compute the future travels distribution.
        The main point is to create the travel distribution according to the
        weight of each produced and attracted travels.
        Takes as input the observed (historical) travels in the form of squared
        matrix,the future produced and attracted ones, respectively, in the
        form of one-line matrices (one for produced and one for attracted)
        and the precision (tolerance).
        Returns a matrix with adjusted travels.
        """

        print()
        print("Enter iter_wgt_dmd method")

        print()
        print("Historical travels matrix travs is, ", travs)
        print()

        # check if the matrices have the correct shape
        # check with the future produced
        if(len(travs) != len(P_is)):
            print("The travels matrix doesn't match with the future\
                  produced! Please fix it.")
            exit()

        # check with the future attracted
        if(len(travs[0]) != len(A_js)):
            print("The travesl matrix doesn't match with the future\
                  attracted! Please fix it.")

            exit()
        
        # function to compare the produced, respectively attracted travels
        # within a certain tolerance
        def comp(s_ih, s_ic, tlr):
            """
            Function within method to compare two values, within tolerance.
            Takes as inputs the lists of to be compared values
            and the precision/tolerance.
            Returns True of False.
            """
            
            # set a flag
            flag = True

            for ih, ic in zip(s_ih, s_ic):
                if(abs(ih - ic) / ih >= tlr): 
                    flag = False
                    break

            return flag
        
        # get produced travels sums on observed travels
        s_Pih = []
        for item in travs:
            s_Pih.append(sum(item))

        # get attracted travels sums on observed travels (cycling on transposes)
        s_Ajh = []   # store the attracted sums

        travs_tt = list(zip(*travs))

        for item in travs_tt:
            s_Ajh.append(sum(item))

        # generate the produced coefficients matrix (produced perspective)
        c_Pi0 = []
        for t, P in zip(travs, s_Pih):
            for trav in t:
                c_Pi0.append(trav/P)
        
        c_Pi = [c_Pi0[i:i + 3] for i in range(0, len(c_Pi0), 3)]

        # generate the attracted coefficients matrix (attracted perspective)
        c_Aj0= []
        for t, A in zip(travs_tt, s_Ajh):
            for trav in t:
                c_Aj0.append(trav/A)

        c_Aj = [c_Aj0[j:j + 3] for j in range(0, len(c_Aj0), 3)]

        print("c_Pi matrix, ", c_Pi)
        print("c_Aj matrix, ", c_Aj)

        cmp_flg = False  # comparison flag to govern the following cycle
        i = 0   # produced passes counter
        j = 0   # attracted passes counter
        
        # allocate initial computed travels matrix
        travsc = []
        for row in travs:
            travsc.append(row)

        while(cmp_flg == False):
                       
            # get produced travels sums on computed travels
            s_Pic = []
            for item in travsc:
                s_Pic.append(sum(item))

            print()
            print("P_is, ", P_is)
            print("s_Pic, ", s_Pic)
            
            cmp_flg = comp(s_Pic, P_is, tlr)
            print(cmp_flg)
            if (comp(s_Pic, P_is, tlr) == False):
                delta_P = []    #list to store the deltas of produced travels
                for Pis, Pic in zip(P_is, s_Pic):
                    delta_P.append(Pis - Pic)
                remind_P = []    #matrix of additions to travels, produced
                for cP, delta_i in zip(c_Pi, delta_P):
                    for c in cP:
                        remind_P.append(c*delta_i)   #append weighted delta
                remind_P = [remind_P[i:i + 3] for i in range(0, len(remind_P), 3)]
                print("remind_P, ", remind_P)
                travsP = []   # list to store adjusted travels matrix, produced
                for remP, trav in zip(remind_P, travsc):
                    for rem, t in zip(remP, trav):
                        travsP.append(rem+t)

                #travsc.clear()

                travsc = [travsP[i:i + 3] for i in range(0, len(travsP), 3)]

                print()
                # print("travs, ", travs)
                # print("travsc, ", travsc)
                
                i += 1

                print()
                print("travsc after pass  i = ", i, "is ", travsc)
            
            s_Pic.clear()

            for item in travsc:
                s_Pic.append(sum(item))

            print()
            print("s_Pic, ", s_Pic)
            
            """
            working on attracted travels
            """

            # transpose de matrix
            travsc_tt = list(zip(*travsc))

            travsc_t = [list(sublist) for sublist in travsc_tt]

            # get attracted travels sums on computed travels (cycling on transposes)
            s_Ajc = []   # store the attracted sums

            for item in travsc_tt:
                s_Ajc.append(sum(item))
            
            print()
            print("A_js, ", A_js)
            print("s_Ajc, ", s_Ajc)

            if (comp(s_Ajc, A_js, tlr) == False):
                delta_A = []    #list to store the deltas of attracted travels
                for Ajs, Ajc in zip(A_js, s_Ajc):
                    delta_A.append(Ajs - Ajc)
                print()
                print("delta_A, ", delta_A)
                print()

                remind_A = []    #matrix of additions to travels, atrracted
                for cA, delta_j in zip(c_Aj, delta_A):
                    for c in cA:
                        remind_A.append(c*delta_j)   #append weighted delta
                remind_A = [remind_A[j:j + 3] for j in range(0, len(remind_A), 3)]
                print("remind_A, ", remind_A)
                print("travsc_tt, ", travsc_tt)
                travsA = []   # list to store adjusted travels matrix, produced
                for remA, trav in zip(remind_A, travsc_tt):
                    #print("remA, ", remA)
                    #print("trav, ", trav)
                    for rem, t in zip(remA, trav):
                        travsA.append(rem+t)

                #travsc.clear()
                #travsc_t.clear()

                travsc_t = [travsA[i:i + 3] for i in range(0, len(travsA), 3)]
                travsc_tt = list(zip(*travsc_t))
                travsc = [list(sublist) for sublist in travsc_tt]
                print("travsA, ", travsA)
                print("travsc_t, ", travsc_t)
                print("travsc_tt, ", travsc_tt)

                print()
                # print("travs, ", travs)
                print("travsc, ", travsc)
                
                j += 1

                print()
                print("travsc_tt after pass j = ", j, "is ", travsc_tt)             
                
            # update the attracted sums
                
            # get attracted travels sums on new computed travels (cycling on transposes)
            s_Ajc.clear()   # clear the computed attracted sums

            for item in travsc_t:
                s_Ajc.append(sum(item))

            print()
            print("s_Ajc, ", s_Ajc)

            # update the produced sums
            s_Pic.clear()

            for item  in travsc:
                s_Pic.append(sum(item))

            print()
            print("s_Pic, ", s_Pic)

            cmp_flgA = comp(s_Ajc, A_js, tlr)
            print("Flag on attracted, ", cmp_flgA)
            
            cmp_flgP = comp(s_Pic, P_is, tlr)
            print("Flag on produced, ", cmp_flgP)

            if(cmp_flgA == True and cmp_flgP == True):
                cmp_flg = True
            else:
                cmp_flg = False

            #if(j == 3):   #block the loop at 3 iterations
            #    cmp_flg = True

            travscr = []     # list to store rounded values, flatten form
            for item in travsc:
                for item in item:
                    travscr.append(round(item))

        #print()
        #print("Final rounded and flatten, ", travscr)
        travscrm = [travscr[i:i + 3] for i in range(0, len(travscr), 3)]
            
        print()
        print("Final rounded matrix, weighted, ", travscrm)
        print("Historical travels matrix, ", travs)
        print("Exit iter_adj_wgt method.")
        
        return travscrm


    # method to compute gravitational model travels projected into the future
    def gravmod_fin(ffs, k_ijs, P_is, A_js):
        """
        Method to compute future travels using gravitational model.
        Takes as inputs the future friction factors matrix, the computed
        calibration coefficients matrix (see next method in this class),
        the matrix of produced travels and the matrix of attracted travels.
        Returns a matrix with future travels determined with gravitational
        model.
        """
        
        # check if the matrices have the same shape
        if(len(k_ijs) != len(ffs) or (len(k_ijs) != len(P_is)) or \
           (len(k_ijs) != len(A_js))):
            print("The matrices doesn't match. Please fix it.")
            exit()
    
        # compute travels with gravitational model
        gvals_fin = []    # to store computed values
        for i in range(len(k_ijs)):
            pdsum = 0
            for j1, j2 in zip(A_js, ffs[i]):
                pdsum = pdsum + j1 * j2
                print("pdsum fin, ", pdsum)
                print("ffs[i] fin, ", ffs[i])
                # print("A_j, ", j1)
            for k1 in range(len(ffs[i])):
                gvals_fin.append((P_is[i] * ffs[i][k1] * A_js[k1] * k_ijs[i][k1] / pdsum))

        print("Future travels obtained with gravitational model, ", gvals_fin)

        # check raw produced travels
        gvals_fin_m0 = [gvals_fin[i:i + 3] for i in range(0, len(gvals_fin), 3)]
        # print(gvals_fin_m0)

    
        # round the number of travels
        gvals_fin_r = []
        for item in gvals_fin:
            gvals_fin_r.append(round(item))
    
        print("Rounded number of furure travels, ", gvals_fin_r)
        print("Total travels sum, ", sum(gvals_fin_r))

        # group flatten list 'gvals_fin_r' as a matrix
        gvals_fin_m = [gvals_fin_r[i:i + 3] for i in range(0,
                         len(gvals_fin_r), 3)]

        print("Matrix of future rounded numbers, ", gvals_fin_m)

        return gvals_fin_r

        
    def ccoeffs(gvalsradj, travs):
        """
        Docstrings goes here ...
        """
        # compute calibration coefficients
        ccoeffs = []
        for row_h, row_c in zip(travs, gvalsradj):
            for t_h, t_c in zip(row_h, row_c):
                ccoeffs.append(round(t_h / t_c, 2))
        ccoeffs_m = [ccoeffs[i:i + 3] for i in range(0, len(ccoeffs), 3)]

        return ccoeffs_m



# *****
"""
Enter the methods calls section here
"""
# gvalsr = GravitMod.gravmod_init(travs, ffs, k_ij0)
# print("gvalsr is, ", gvalsr)

# gvalsr_m = [gvalsr[i:i + 3] for i in range(0, len(gvalsr), 3)]

#gvalsadjA = GravitMod.iter_adj_in(travs, gvalsr)
#gvalsadjB = GravitMod.iter_adj_wgt(travs, gvalsr)

#ccoeffsA = GravitMod.ccoeffs(gvalsadjA, travs)
#ccoeffsB = GravitMod.ccoeffs(gvalsadjB, travs)

travsc_furn = GravitMod.furness(travs, P_is, A_js)
travsc_wgtd = GravitMod.iter_wgt_dmd(travs, P_is, A_js)
travs_frat = GravitMod.fratar(travs, P_is, A_js)
print()
print("Matrix of travels obtained with weighted coefficients is, ",
      travsc_wgtd)

#print("Adjusted matrix A, ", gvalsadjA)
#print("Adjusted matrix B, ", gvalsadjB)
#print("Calibration coefficients matrix A, ", ccoeffsA)
#print("Calibration coefficients matrix B, ", ccoeffsB)

# print("ccoeffs, ", ccoeffs)

# ccoeffs_it = GravitMod.ccoeffs(gvalsradj_it, travs)

# print("ccoeffs_it, ", ccoeffs_it)


# print("Calibration coefficients, ", ccoeffs_m)

# gvalsr_fin = GravitMod.gravmod_fin(ffs_f, ccoeffs_m, P_is, A_js)

# print("Future number of rounded travels, ", gvalsr_fin)

# u_a, u_t = modopt(tca, tct, tda, tdt)

# w_a, w_t = logit(u_a, u_t)
