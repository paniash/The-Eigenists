def QPE(thet, sz):
    qc = QuantumCircuit(5)
    qc.x(0)
    qc.u1(1, 0)
    qc.x(0)
    qc.u1(1, 0)
    qc.u3(-2*thet, -np.pi/2, -np.pi/2, 0)
    
    # Gates for A3
    for i in range(0, sz-1):
        q_controls = []
        qc.cx(i, i+1)
        q_controls.append(i+1)
        for j in range(i, 0, -1):
            qc.cx(i, j-1)
            q_controls.append(j-1)
    #         Multicontrolled x rotation
        if(len(q_controls)>1):
            qc.mcmt(CU3Gate(-2*thet, -np.pi/2, -np.pi/2), q_controls, i)
        else:
            qc.cu3(-2*thet, -np.pi/2, -np.pi/2, q_controls[0], i)
    #         # Uncompute
        for j in range(0, i):
            qc.cx(i, j)
        qc.cx(i, i+1)
    
    return qc