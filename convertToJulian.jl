function convertToJulian(D, M, Y, h, m, s)
JD = 367*Y - floor(7*(Y+floor((M+9)/12))/4) + floor(275*M/9) + D + 1721013.5 + ((s/60 + m)/60 + h)/24  
MJD = JD - 2400000.5
end

