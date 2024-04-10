function [ DFP ] = findMatrixFPrime( mu, p )
    global DRXA XAP DRXB XBP FP0;
    DFP = zeros(0);
    
    DFP = DRXA*XAP+DRXB*XBP;
    DFP = DFP^(-1);
    DFP = -DFP;
    DFP = DFP*FP0;
    DFP = double(DFP);
end
