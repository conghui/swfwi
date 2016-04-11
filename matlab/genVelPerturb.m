function genVelPerturb(NXMAX, NZMAX, numVel, maxPerturb)

if ischar(NXMAX)
    NXMAX = str2num(NXMAX);
    NZMAX = str2num(NZMAX);
    numVel = str2num(numVel);
    maxPerturb = str2num(maxPerturb);
end

filename = ['vp' '-nx' num2str(NXMAX) '-nz' num2str(NZMAX) '-nv' num2str(numVel) '-p' num2str(maxPerturb) '.bin']
fid = fopen(filename, 'wb');

NFre=20;
NXMAX = NXMAX + NFre;
NZMAX = NZMAX + NFre;

rng('default');
rng(1);             % Set the seed

for i = 1:numVel
    matSpec=zeros(NXMAX,NZMAX);
    for indX=1:NXMAX
        for indZ=1:NZMAX
            if ((indX<NXMAX-NFre && indX>NFre*NXMAX/NZMAX)||(indZ<NZMAX-NFre && indZ>NFre))
                continue;
            end
            matSpec(indX,indZ)=(2*rand(1)-1)+(2*rand(1)-1)*1i;
        end
    end

    matSpace = real(ifft2(matSpec));

    % get the actual size without boundary
    perturb = matSpace(11:end - 10, 11:end - 10);

    % find out the maxsimum absolute value
    maxabs = max(max(abs(perturb)));

    velper = maxPerturb / maxabs * perturb;

    fwrite(fid, velper', 'float32');

    %mesh(velper);
end % end for

fclose(fid);

end % end function
