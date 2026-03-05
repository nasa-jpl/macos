
% Compare OPDs for ZernGridData surface testing

% Zernike comparison

load -ascii tst_ZGDzern2_ZGD.txt
opdZzgd = tst_ZGDzern2_ZGD;

load -ascii tst_ZGDzern2_zern.txt
opdZzern = tst_ZGDzern2_zern;

dzern = opdZzgd-opdZzern;
nnzz = nnz(dzern)

figure(1)
clf
set(1,'pos',[680   678   804   320])
dimage([opdZzgd opdZzern dzern])
title('Comparison of ZrnGrData and Zernike Surface Types')
xlabel(['nnz(opdZzgd-opdZzern) = ' num2str(nnzz)])


% Data comparison

load -ascii tst_ZGDdata2_ZGD.txt
opdDzgd = tst_ZGDdata2_ZGD;

load -ascii tst_ZGDdata2_grid.txt
opdDdata = tst_ZGDdata2_grid;

ddata = opdDzgd-opdDdata;
nnzd = nnz(ddata)

figure(2)
clf
set(2,'pos',[619   576   804   320])
dimage([opdDzgd opdDdata ddata])
title('Comparison of ZrnGrData and GridData Surface Types')
xlabel(['nnz(opdZzgd-opdZzern) = ' num2str(nnzd)])


