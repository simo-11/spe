% reproduce bug mentioned in support case 07309712
fn='echodemo_bug.m';
fd = fopen(fn, 'w');
for sn=1:2
    fprintf(fd,'%%%% %d.\nfprintf("%d. section");\n',sn,sn);
    fclose(fd);
    echodemo('echodemo_bug',sn);
    fd=fopen(fn,"a");
end