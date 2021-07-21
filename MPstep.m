% simple 2D step finding algorithm for a fixed interaction potential
% trajectory by Marek Piliarik and Lukasz Bujak for manuscript submitted to
% Small Methods "Fast Leaps Between Millisecond Confinements Govern Ase1 
% Diffusion Along Microtubules"
% Sample data available from Figshare - link will be available upon
% aceptance.

%% load the data
Ase=load ('all57.mat') ;

    T=Ase.t';
    X=Ase.x'/8; % in tubulin dimers
    Y=Ase.yprime'/6; % in tubulin dimers
 

    
 %% smooth and make histogram
    Cpos=X+1i*Y;
    SM=15;
    SM2=SM;
    MinDisc=0.0;
    MinDwell=round(44*10^MinDisc/2);
    sCpos=smooth(Cpos,SM);
  figure(2);
  discr=0.05;
%  EdgeX=-40:discr:40;
  EdgeX=-90:discr:60;
  EdgeY=-15:discr:15;
  hhh=histogram2(real(sCpos),imag(sCpos),EdgeX,EdgeY);
  hhh=imsmooth(hhh.Values,0.2/discr);
  fmap((hhh));hold on; plot (smooth(Cpos-EdgeX(1)-EdgeY(1)*1i+discr*0.5*(1+1i),SM)/discr);hold off
    
  
%% index peaks
histsurf=hhh;
cdata=smooth(Cpos,SM2);
template=histsurf.*0;
template(1:end-1,:)=histsurf(1:end-1,:)>histsurf(2:end,:);
template(2:end,:)=template(2:end,:) & histsurf(2:end,:)>histsurf(1:end-1,:);
template(:,1:end-1)=template(:,1:end-1) & histsurf(:,1:end-1)>histsurf(:,2:end);
template(:,2:end)=template(:,2:end) & histsurf(:,2:end)>histsurf(:,1:end-1);
template(template(:)>0)=1:sum(template(:));
[~, I]=sort(histsurf(:),'descend');

for ii=1:length(I)
    if histsurf(I(ii))==0
        break
    end
     [xx, yy]=ind2sub(size(histsurf),I(ii));
     if template(I(ii))==0
         neighbors=[];
         nc=1;
        for xxx=xx-1:xx+1
            for yyy=yy-1:yy+1
         if (xxx>0)&&(yyy>0)&&(xxx<size(template,1))&&(yyy<size(template,2))&& ...
                 template(xxx,yyy)>0 %histsurf(xx,yy)< histsurf(xxx,yyy);
             template(I(ii))=template(xxx,yyy);
             neighbors(nc)=template(xxx,yyy);
             nc=nc+1;
         end
            end
        end
        if max(neighbors)>min(neighbors)
            template(I(ii))=0;
        end
     end
     
end
  template(:,2:end-1)=round(sqrt(template(:,2:end-1).*sqrt(template(:,3:end).*template(:,1:end-2))));
  template(2:end-1,:)=round(sqrt(template(2:end-1,:).*sqrt(template(3:end,:).*template(1:end-2,:))));
%}
%h1=fspecial('gauss',9,4);h1(:)=h1(:)>max(h1(:))/2;
%template=imfilter(template,h1);

%% plot segmentation
  figure(3);
  subplot(2,2,1);
  
fmap(template+1); colormap colorcube
;hold on; plot (smooth(Cpos-EdgeX(1)-EdgeY(1)*1i,50)/discr,'w');hold off
idata=zeros(length(cdata),1);
for ii=1:length(cdata)
    EdgeXi=find(diff( EdgeX>real(cdata(ii)) ));
    EdgeYi=find(diff( EdgeY>imag(cdata(ii)) ));
    if isempty(EdgeYi)
        EdgeYi=length(EdgeY)-1;
    end
    idata(ii)=template(EdgeXi,EdgeYi);
end


for ii=2:length(idata)
    if idata(ii)==0
        idata(ii)=idata(ii-1);
    end
    if ii>1+MinDwell
%         if idata(ii+MinDwell)==idata(ii)
%             idata(ii:ii+MinDwell)=idata(ii);
%         end
        if sum(abs(diff(idata(ii-MinDwell:ii-1))))>1
            SINDX=round( ...
                sum( (1:MinDwell)' .* (abs(diff(idata(ii-MinDwell:ii)))) ) ...
                / sum( abs(diff(idata(ii-MinDwell:ii))) ) ...
                        );
            idata(ii-MinDwell:ii-MinDwell+SINDX)=idata(ii-MinDwell);
            idata(ii-MinDwell+SINDX:ii)=idata(ii);
        end
    end
end


figure(3);
subplot(2,2,2);
scatter(real(cdata),imag(cdata),3,idata);
colormap colorcube
%idata=cumsum(abs(diff(idata)));
cstepdata=cdata*0;
for ii=1:max(idata)
    cstepdata(idata==ii)=mean(cdata(idata==ii));
end
  figure(3);
  subplot(2,2,[3 4]);
plot(T,[X, Y]);hold on;
plot(T,[real(cstepdata), imag(cstepdata)],'LineWidth',3);hold off
%title(sprintf('trace #%d',jj));
% saveas(gcf,sprintf('Ase%d_Confinements.tif',jj));
% saveas(gcf,sprintf('Confinements_Ase%d.fig',jj));



%% displacement and dwell-time histogram
    stepsX=diff(cstepdata(:));
    stepsTX=T(abs(stepsX)>0);
    stepsX=(cstepdata(1+1:end)-cstepdata(1:end-1));
    stepsX=stepsX(abs(stepsX)>0);
    hh=histcounts2(real(stepsX),imag(stepsX),-5:.05:5,-5:.05:5);
    HHX=hh>0;
  figure(4);
  subplot(1,2,1);

    fmap(imsmooth(double(HHX),5));axis equal; colorbar off
  figure(4);
  subplot(1,2,2);
    hhbin=10.^(MinDisc:0.1:3);
    hh=histogram(diff(stepsTX),hhbin,'Normalization','pdf');
    HHTX=hh.Values;

    axis([10^MinDisc 1000 0 max((HHTX))])
    set(gca,'YScale','lin')
    set(gca,'XScale','log')
% saveas(gcf,sprintf('Ase%d-stepstat.tif',jj));
% saveas(gcf,sprintf('stepstat-Ase%d.fig',jj));


function fmapres = fmap(datain, zmin, zmax)

if nargin<3
  zmin=0;
  zmax=0;
end;

if nargin==2
  figure(zmin)
end;
data(:,:)=datain;
pcolor(double(data)'); 
caxis([zmin zmax]);
shading flat; colormap hot; colorbar;
cmapres=size(data);
end

function ims=imsmooth(data, width, verbose)
if (nargin<2)
  width=2;
end;
if (nargin<3)
  verbose=true;
end;
clear h1;
h1=fspecial('gaussian',ceil(6*width), width);

if verbose 
    fprintf('smoothing images ------');
end
ims=data;
for ii=1:length(data(1,1,:))
  ims(:,:,ii)=imfilter(ims(:,:,ii),h1,'replicate','conv');
  if (floor(ii/10)==ii/10) && verbose  
    fprintf('\b\b\b\b\b\b\b %06d', ii);
  end;
end;
% fprintf('\n');
end