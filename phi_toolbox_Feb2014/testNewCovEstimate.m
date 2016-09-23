
%instead of generating data will use the eeg data

load('A:\drorcNew\Current\Data\FrancescaSleepData\data_chansKept_centeralANDNOTFlat_aveReref_detrended_fltrd_0.5_35\Subject2\N1\data.mat')

%% first just calculate the covs

maxChans=[40 50 60 70];

for maxChanIndx=1:length(maxChans)
    maxChan=maxChans(maxChanIndx);
    tic;
    
    for epochIndx=1:length(data)

    %%% just to see some performace issues
         tOrig=tic ;
        [CovXt{maxChanIndx,epochIndx} CovXtXtau{maxChanIndx,epochIndx} CovXtau{maxChanIndx,epochIndx}] = Cov_comp_sample(data(epochIndx).data(1:maxChan,:),10);
         tOrigEpoch(maxChanIndx,epochIndx)=toc(tOrig);
    %     disp(['cond cov sample ' num2str(cond(CovXt))]);

         tShrink=tic; 
        [Cov_X2{maxChanIndx,epochIndx} Cov_XY{maxChanIndx,epochIndx} Cov_YX{maxChanIndx,epochIndx} Cov_Y2{maxChanIndx,epochIndx}] = Cov_comp_shrink(data(epochIndx).data(1:maxChan,:),10,0);
        tShrinkEpoch(maxChanIndx,epochIndx)=toc(tShrink);


    end

    %% 
    
    

    for epochIndx=1:length(data)

        condNumOrig(maxChanIndx,epochIndx) = cond(CovXt{maxChanIndx,epochIndx});
        condNumShrink(maxChanIndx,epochIndx) = cond(Cov_X2{maxChanIndx,epochIndx});

    end



    %% and calc phi
    for epochIndx=1:length(data)

        [phi_star{maxChanIndx,epochIndx} phi_star_fixed{maxChanIndx,epochIndx} phi_MI{maxChanIndx,epochIndx} phi_MI_fixed{maxChanIndx,epochIndx} MI{maxChanIndx,epochIndx} MI_fixed{maxChanIndx,epochIndx} phi_H{maxChanIndx,epochIndx} beta{epochIndx} betaFixed] = ...
            phi_comp(CovXt{maxChanIndx,epochIndx},CovXtXtau{maxChanIndx,epochIndx},CovXtau{maxChanIndx,epochIndx},eye(maxChan),1:maxChan,1);

        %also check masafui's new version  
        disp 'new'
        [phi_star2{maxChanIndx,epochIndx} phi_star_fixed2{maxChanIndx,epochIndx} phi_MI2{maxChanIndx,epochIndx} phi_MI_fixed2{maxChanIndx,epochIndx} MI2{maxChanIndx,epochIndx} MI_fixed2{maxChanIndx,epochIndx} phi_H2{maxChanIndx,epochIndx} beta{maxChanIndx,epochIndx} betaFixed] = ...
            phi_comp(Cov_X2{maxChanIndx,epochIndx},Cov_XY{maxChanIndx,epochIndx},Cov_Y2{maxChanIndx,epochIndx},eye(maxChan),1:maxChan,1);

        disp(epochIndx);
    end
    
end

%% now plot run time

figure
plot(maxChans,mean(tOrigEpoch,2),'b+-',maxChans,mean(tShrinkEpoch,2),'r+-')
title('run time comparison')
ylabel('time sec');
xlabel('num chans');
legend({'sample cov', 'shrink cov'});

figure
plot(maxChans,mean(condNumOrig,2),'b+-',maxChans,mean(condNumShrink,2),'r+-')
title('condition number obtained using cond()')
ylabel('cond number (ratio of largest to smallest singular vals');
xlabel('num chans');
legend({'sample cov', 'shrink cov'});

figure
array_phi_stars_orig=vertcat(phi_star{2,:});
array_phi_stars_shrink=vertcat(phi_star2{2,:});


array_phi_starsFixed_orig=vertcat(phi_star_fixed{2,:});
array_phi_starsFixed_shrink=vertcat(phi_star_fixed2{2,:});
    

subplot(2,2,1)
plot(1:length(array_phi_stars_orig),array_phi_stars_orig(:,1),1:length(array_phi_stars_orig),array_phi_stars_shrink(:,1),'r')
title('phi')
xlabel('trial');
legend({'phi star with sample cov', 'phi star with shrinkage'}) 

subplot(2,2,3)
plot(1:length(array_phi_stars_orig),array_phi_stars_orig(:,2),1:length(array_phi_stars_orig),array_phi_stars_shrink(:,2),'r')
title('phi NORMALIZED')
xlabel('trial');
%%%%% phistar fixed

subplot(2,2,2)
plot(1:length(array_phi_stars_orig),array_phi_starsFixed_orig(:,1),1:length(array_phi_stars_orig),array_phi_starsFixed_shrink(:,1),'r')
xlabel('trial');
title('phi fixed')


subplot(2,2,4)
plot(1:length(array_phi_stars_orig),array_phi_starsFixed_orig(:,2),1:length(array_phi_stars_orig),array_phi_starsFixed_shrink(:,2),'r')
xlabel('trial');
title('phi fixed NORMALIZED')

    
suptitle('phis calculated with 60 chans');

%% the phi star fixed spikes occur at epoch 12 18 and 22

figure
subplot(1,3,1);
plot(data(12).data(1:maxChan,:)')

subplot(1,3,2);
plot(data(18).data(1:maxChan,:)')

subplot(1,3,3);
plot(data(22).data(1:maxChan,:)')

suptitle('data for trials with spikes in Integrated Information calc');  


