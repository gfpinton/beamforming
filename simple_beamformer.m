%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GIANMARCO PINTON
% WRITTEN: 2020-09-19
% LAST MODIFIED: 2020-09-19
% BASIC TABLE-BASED BEAMFORMING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load rf_pointtarget.mat % pxducer: transducer data (single precision)
			% dX: grid spacing in pxducer (m)
			% dT: time samplilng in pxducer (s)
			% txD: time delays for transmit profiles (time pixels)
%%% Basic variable definition %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c0=1540;
%%% Image grid defintion %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
deps = 10e-3:0.125e-3/4:0.9*size(pxducer,1)*dT*c0/2;
lats = -5e-3:0.25e-3/4:5e-3;
xducercoords=[(0:size(pxducer,2)-1)' zeros(size(pxducer,2),1)];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
bm=zeros(length(lats),length(deps),'single');
idps=zeros(length(lats),length(deps),size(pxducer,2)+1,'int64');

idt0=max(txD);

for ii=1:length(lats)
  lat=lats(ii);
  for jj=1:length(deps)
    dep=deps(jj);
    fcen=([lat/dX+mean(xducercoords(:,1)) dep/dX ]);
    [txdd txmdd]=focusProfile2(fcen,mean(xducercoords),dT/dX*c0);
    [rxdd rxmdd]=focusProfile2(fcen,xducercoords,dT/dX*c0);
    idx=find(xducercoords(:,1)>-inf); % no fnumber restriction here
	   % the total time delay is the initial offset (idt0) +
	   % the transmit delay (txmdd) +
	   % the receive profile (rxmdd) 
    idp=int64(round(((size(pxducer,1)*(idx-1))+idt0+txmdd+rxmdd)));
    idps_mat(ii,jj,1)=length(idp);
    idps_mat(ii,jj,2:length(idp)+1)=idp;
  end
end
%% out of bounds check %
idm=find(idps_mat<1 | idps_mat>size(pxducer,1)*size(pxducer,2));
if(idm)
  length(idm)/length(idps_mat(:))
  idps_mat(idm)=1;
end

tic
for ii=1:length(lats)
  for jj=1:length(deps)
    bm(ii,jj)=sum(pxducer(idps_mat(ii,jj,2:idps_mat(ii,jj,1)+1)));
  end
end
toc

imagesc(db(abs(hilbert(bm'))))

pxducer2=pxducer; pxducer2(:,:,2)=pxducer;
tic
bmf=mex_beamform_omp4_float(int64(idps_mat),single(pxducer2));
toc
imagesc(db(abs(hilbert(squeeze(bmf(:,1,:))))))
