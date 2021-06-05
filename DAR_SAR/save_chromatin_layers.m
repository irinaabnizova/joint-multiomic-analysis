function [add]=save_chromatin_layers(DAR,outFilename)
%==========write DAR or SAR with GC levels for each lineage
%------------------------Input
 %SAR=[chr st en chp' lev_ES' lev_EnS' lev_CS' ind'];
    %    1    2  3  4     5        6     7       8
    add=25;
    chr=DAR(:,1);
    starts=DAR(:,2)-add;
    ends=DAR(:,3)+add;
    ch=DAR(:,4);
    lev_EcS=DAR(:,5);
    lev_EnS=DAR(:,6);
    lev_MS=DAR(:,7);
    ind=DAR(:,8);
    
       
  %--------------------DARs/SARs 8 fileds
  %textFilename=sprintf('SARL_EctEndMes_c25_aL02_GCPu.txt');
  fp = fopen(outFilename, 'w');
  for i=1:length(starts),
       fprintf(fp,'%d\t%d\t%d\t%.2f\t%.2f\t%.2f\t%.2f\t%d\n',chr(i),starts(i),ends(i),ch(i),lev_EcS(i),lev_EnS(i),lev_MS(i),ind(i));
  end
  fclose(fp);
 
  