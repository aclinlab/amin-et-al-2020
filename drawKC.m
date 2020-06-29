i = 1847;

kc = neurSkel(strcat('KCskel',num2str(i),'.csv'));
kc = kc.getSynSets({strcat('KCoutputs',num2str(i),'.csv'), ...
    strcat('KCinputs',num2str(i),'.csv')}, [1 2]);
kc = kc.mapOrigSynsToSynLocs;
kc.drawKCSkeleton();
