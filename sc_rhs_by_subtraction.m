%% settings
%ao.models=["poly44","cubicinterp","tps"];
ao.models=["cubicinterp","linearinterp"]; %#ok<*NBRAK2>
%ao.cubs=["integral2","glaubitz","rbfcub"];
ao.cubs=["integral2"];
ao.debugLevel=0;
ao.plot=0;
ao.scat_type='halton';
ao.cards=[50,100,400,500];
ao.rsquareMin=0.9;
%number of nodes in section properties model
ao.n=["122","682","1319","1453";
      "86", "460","848", "995"];
%ao.n="375";
ao.oms=[4,10,20,100];% output mesh size, number of points in each direction
add_lib_to_path
E=1.e9;
nu=0.3;
G=E/(2*(1+nu));
L=2;
wa=[160;140];% mm
ha=[80;60];% mm
%% interpolate based on section properties results
ns=size(ao.n,2);
r=cell(2*ns,1);
for j=1:2
    h=ha(j);
    w=wa(j);
    for i=1:ns
        expected_name=sprintf("gen/warping-rectangle-%g-%g-%s.csv",...
            h,w,ao.n(j,i));
        if ~isfile(expected_name)
            error('results are not available, expected to have file %s',...
                expected_name);
        end
        r{(j-1)*ns+i}=testRectangle(height=h,width=w,models=ao.models,...
        cubs=ao.cubs,debug=ao.debugLevel,cards=ao.cards,n=ao.n(j,i),...
        latex=0);
    end
end    
%% write warping results based on interpolation f
ns=size(ao.n,2);
for j=1:2
    h=ha(j);
    w=wa(j);
    for i=1:ns
        fr=r{(j-1)*ns+i}{1};
        list=dir(fr.file);
        expected_name=sprintf("warping-rectangle-%g-%g-%s.csv",...
            h,w,ao.n(j,i));
        if list.name~=expected_name
            error('results are for %s,\nbut current model expects %s',...
                list.name,expected_name);
        end
        ms=size(ao.models,2);
        for mi=1:ms
            model=ao.models(mi);
            es=sprintf("f=fr.%s_fit;",model);
            eval(es);
            omss=size(ao.oms,2);
            for omi=1:omss
                oms=ao.oms(omi);
                es=sprintf("fn='%s-%g-%g-%s-%dx%d.xlsx';",...
                    model,h,w,ao.n(j,i),oms,oms);
                eval(es);
                save_interpolation_results(f,w/1000,h/1000,oms,fn);
            end
        end
    end
end
