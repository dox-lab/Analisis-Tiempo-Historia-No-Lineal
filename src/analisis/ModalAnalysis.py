def ModalAnalysis(numEigen, pflag=1, outname=None):
    """
    This script returns the modal properties of an OpenSeespy model.

    Parameters
    ----------
    numEigen : int
        Number of eigenvalues to calculate.
    pflag    : int (1 or 0)
        Flag to print output information on screen.
    outname  : str, optional
        If provided and pflag == 1, the modal properties will be saved to outname.csv.

    Returns
    -------
    T        : numpy.ndarray
        Period array for the first numEigen modes.
    Mratios  : dict
        Effective modal mass participation ratios.
    Mfactors : dict
        Modal participation factors.
    Mtots    : dict
        Total activated masses.
    """
    import numpy as np
    import openseespy.opensees as ops
    import sys

    ops.wipeAnalysis()
    ops.numberer("Plain")
    ops.system('FullGeneral')
    ops.algorithm('Linear')
    ops.analysis('Transient')

    # Extract the mass matrix (unrestrained part)
    ops.integrator('GimmeMCK', 1.0, 0.0, 0.0)
    ops.analyze(1, 0.0) 
    N = ops.systemSize()
    Mmatrix = ops.printA('-ret')
    Mmatrix = np.array(Mmatrix)
    Mmatrix.shape = (N, N)
    print('\nExtracting the mass matrix, ignore the warnings...')
        
    # Determine maximum number of DOFs per node
    NDF = 0
    for node in ops.getNodeTags():
        temp = len(ops.nodeDOFs(node))
        if temp > NDF:
            NDF = temp

    DOFs = []
    used = {}
    ldict = {}
    Mratios = {}
    Mfactors = {}
    
    for i in range(1, NDF + 1):
        ldict[i] = np.zeros([N, 1])
        Mratios[i] = np.zeros(numEigen)
        Mfactors[i] = np.zeros(numEigen)

    idx = 0
    for node in ops.getNodeTags():
        used[node] = []
        ndof = len(ops.nodeDOFs(node))
        for j in range(ndof):
            temp = ops.nodeDOFs(node)[j]
            if temp not in DOFs and temp >= 0:
                DOFs.append(temp)
                used[node].append(j + 1)
                ldict[j + 1][idx, 0] = 1
                idx += 1

    Mmatrix = Mmatrix[DOFs, :][:, DOFs]

    # Calculate the total mass for each direction
    Mtots = {i: (ldict[i].T @ Mmatrix @ ldict[i])[0, 0] for i in range(1, NDF + 1)}

    # Perform eigenvalue analysis
    ops.wipeAnalysis()
    listSolvers = ['-genBandArpack', '-fullGenLapack', '-symmBandLapack']
    ok = 1
    for s in listSolvers:
        print(f"Using {s[1:]} as solver...")
        try:
            eigenValues = ops.eigen(s, numEigen)
            ok = 0
            break
        except:
            pass

    if ok != 0:
        print("Error in eigenvalue analysis.")
        sys.exit()
    else:
        Lambda = np.asarray(eigenValues)
        Omega = np.sqrt(Lambda)
        T = 2 * np.pi / Omega
        frq = 1 / T

    for mode in range(1, numEigen + 1):
        idx = 0
        phi = np.zeros([N, 1])
        for node in used:
            for dof in used[node]:
                phi[idx, 0] = ops.nodeEigenvector(node, mode, dof)
                idx += 1

        phi /= np.sqrt(phi.T @ Mmatrix @ phi)
        Mn = phi.T @ Mmatrix @ phi

        for j in range(1, NDF + 1):
            if Mtots[j] != 0:
                Ln = phi.T @ Mmatrix @ ldict[j]
                Mnstar = (Ln ** 2 / Mn)[0, 0]
                Mfactors[j][mode - 1] = Ln / Mn
                Mratios[j][mode - 1] = Mnstar / Mtots[j] * 100

    for j in range(1, 7):
        if j not in Mratios:
            Mratios[j] = np.zeros(numEigen)
            Mfactors[j] = np.zeros(numEigen)

    del Mratios[6], Mratios[5], Mratios[4]
    del Mfactors[6], Mfactors[5], Mfactors[4]

    sM1 = np.cumsum(Mratios[1])
    sM2 = np.cumsum(Mratios[2])
    sM3 = np.cumsum(Mratios[3])

    if pflag == 1:
        arguments = []
        arguments.append('Modal Periods and Frequencies')
        arguments.append('%4s|%8s|%10s|%12s|%12s' % ('Mode', 'T [sec]', 'f [Hz]', '\u03C9 [rad/sec]', '\u03BB [rad\u00b2/sec\u00b2]'))
        for mode in range(numEigen):
            arguments.append('%4s|%8s|%10s|%12s|%12s' % (
                "{:.0f}".format(mode + 1), "{:.4f}".format(T[mode]), "{:.3f}".format(frq[mode]), 
                "{:.2f}".format(Omega[mode]), "{:.2f}".format(Lambda[mode])))

        arguments.append('Total Activated Masses')
        arguments.append('%8s|%8s|%8s' % ('M\u2081', 'M\u2082', 'M\u2083'))
        arguments.append('%8s|%8s|%8s' % ("{:.2f}".format(Mtots[1]), "{:.2f}".format(Mtots[2]), "{:.2f}".format(Mtots[3])))
        
        arguments.append('Modal Mass Participation Factors')
        arguments.append('%4s|%7s|%7s|%7s' % ('Mode', '\u0393\u2081', '\u0393\u2082', '\u0393\u2083'))
        for mode in range(numEigen):
            arguments.append('%4s|%7s|%7s|%7s' % (
                "{:.0f}".format(mode + 1), "{:.3f}".format(Mfactors[1][mode]), 
                "{:.3f}".format(Mfactors[2][mode]), "{:.3f}".format(Mfactors[3][mode])))
        
        arguments.append('Effective Modal Mass Participation Ratios [%]')
        arguments.append('%4s|%7s|%7s|%7s' % ('Mode', 'U\u2081', 'U\u2082', 'U\u2083'))
        for mode in range(numEigen):
            arguments.append('%4s|%7s|%7s|%7s' % (
                "{:.0f}".format(mode + 1), "{:.3f}".format(Mratios[1][mode]), 
                "{:.3f}".format(Mratios[2][mode]), "{:.3f}".format(Mratios[3][mode])))
        
        arguments.append('Cumulative Effective Modal Mass Participation Ratios [%]')
        arguments.append('%4s|%7s|%7s|%7s' % ('Mode', '\u2211U\u2081', '\u2211U\u2082', '\u2211U\u2083'))
        for mode in range(numEigen):
            arguments.append('%4s|%7s|%7s|%7s' % (
                "{:.0f}".format(mode + 1), "{:.3f}".format(sM1[mode]), 
                "{:.3f}".format(sM2[mode]), "{:.3f}".format(sM3[mode])))

        arguments = '\n'.join(arguments)
        print(arguments)

        if outname is not None:
            with open(outname + '.csv', 'w', encoding='utf-32') as f:
                f.write(arguments)

    return T, Mratios, Mfactors, Mtots
