% calc_frap - makes the actual analysis of the FRAP data
%   [k,D1k,D1,D2,gamma2,gamma0,v_conv,a_conv]=...
%   calc_frap(file_in,dx,Ipre,Idark,istart,t,trc,Rp_max,nr) where:
%
%   k =  spatial frequencies [um^-1]
%   D1k = D1(k2), the diffusion coefficient determined for each k
%       separately [um^2s^-1]
%   D1 and D2 = effective values of the diffusion coefficients D1 and
%       D2 [um^2s^-1]
%   gamma2 = the intensity fraction of component 2
%   gamma0 = the intensity fraction of immobile molecules
%   v_conv = the absolute velocity of the tracked center of mass [um/s]
%   a_conv = the angle the tracked center of mass moves in, with 0 degrees
%       corresponding to moving horizontally from the left to the right
%       [degrees]
%
%   file_in = the filename + directory of the image stack
%   dx = the pixel size (square pixels assumed)
%   Ipre = the pre-bleach intensity
%   Idark = the dark count intensity
%   istart = the first post-bleach frame
%   t = the times for each frame
%   trc = parameter that if equal to 'y' leads to tracking of the center of
%       mass for each frame (otherwise the center of mass is determined from
%       the first post-bleach frame)
%   Rp_max = the maximum radial value used in the analysis [pixels]
%   nr = a vector for all frames in the image stack. If nr==1 then this
%       image will be part of the analysis.

function [k,D1k,D1,D2,gamma2,gamma0,v_conv,a_conv]=...
    calc_frap(file_in,dx,Ipre,Idark,istart,t,trc,Rp_max,nr)

% closes the figures
close(figure(1))
close(figure(2))

% Meshpoints
x=(1:1:size(Ipre,2));
y=(size(Ipre,1):-1:1);
[Xx,Yy]=meshgrid(x,y);

% Calculates the centre of mass from the first frame after bleaching
I=double(imread(file_in,istart))-Idark;
I2=I./Ipre;

fig=figure(1);
set(fig,'Color','w');
I2_min=max(0.5,min(abs(I2(:))));
[I2_max,max_pos]=max((I(:)+Idark)./(Ipre(:)+Idark));
I2ve=I2(:);
I2_max=I2ve(max_pos);

warning off;
imshow(I2,[I2_min,I2_max]);
ax1 = gca;
set(ax1,'FontSize',8)
axis square
h=title({'Define a polygon around the bleached spot (close with the right mouse button).';...
    'Exit by double clicking the inner part of the polygon.'},...
    'Color','k','BackgroundColor','w');
hp=get(h);
pos=hp.Position;
pos(2)=pos(2)+50;
set(h,'Position',pos);
BW=roipoly; % the region in which the bleached spot is situated

Xv=Xx(BW);
Yv=Yy(BW);
Iv=I2(BW)/mean(I2(:));
 
% Calculates the centre of mass (x_cm,y_cm)
options=optimset('Display','off');

pcm0=[1,1-min(Iv),mean(Xv),mean(Yv),(max(Xv)-min(Xv))/4];
pcm=lsqnonlin(@(p) gauss_fkn(p,Xv,Yv,Iv),pcm0,...
    [0,-Inf,min(Xv),min(Yv),0],[Inf,Inf,max(Xv),max(Yv),Inf],options);

x_cm=pcm(3);
y_cm=pcm(4);
pcm0=pcm;

%x_cm=sum(sum(Xx.*(1-I2).*BW))/sum(sum((1-I2).*BW));
%y_cm=sum(sum(Yy.*(1-I2).*BW))/sum(sum((1-I2).*BW));
Rr=sqrt((Xx-x_cm).^2+(Yy-y_cm).^2);
wid=max(Rr(BW==1)); % max r-value to use to calculate the center of mass

if ~isempty(Rp_max)
    BWframe=zeros(size(I));
    BWframe(Xx>(x_cm-Rp_max) & Xx<(x_cm+Rp_max) & Yy>(y_cm-Rp_max) & Yy<(y_cm+Rp_max))=1;
    x=x(x>(x_cm-Rp_max) & x<(x_cm+Rp_max));
    y=y(y>(y_cm-Rp_max) & y<(y_cm+Rp_max));
else
    BWframe=ones(size(I));
end

x_cm=x_cm-x(1)+1;
x=x-x(1)+1;
y_cm=y_cm-y(end)+1;
y=y-y(end)+1;
[Xx,Yy]=meshgrid(x,y);

Ipre=reshape(Ipre(BWframe==1),[length(y),length(x)]);
I2=reshape(I2(BWframe==1),[length(y),length(x)]);
BW=reshape(BW(BWframe==1),[length(y),length(x)]);

figure(1)
imshow(I2,[I2_min,I2_max]);
hold on
plot3([x_cm-10,x_cm+10],[max(Yy(:))-y_cm+1,max(Yy(:))-y_cm+1],max(I2(:))*ones(1,2),...
    [x_cm,x_cm],[max(Yy(:))-y_cm+1-10,max(Yy(:))-y_cm+1+10],max(I2(:))*ones(1,2),'LineWidth',2,'Color','red')
hold off
pause(0.2)

% Defines the options for nonlinear fitting
options=optimset('Display','off');

Nrp=512;    % the length of r

% Determines the radial distance r and the radius of the field of view R
[test,r,R,Rp_max]=I_radial(Nrp,x_cm,y_cm,Xx,Yy,Xx,Rp_max);
r=r*dx;
R=R*dx;

% Initiates variables
wf=zeros(sum(nr==1),1);
x_cm2=zeros(sum(nr==1),1);
y_cm2=zeros(sum(nr==1),1);
Ixx=zeros(size(x_cm2));
Iyy=zeros(size(x_cm2));
Ixy=zeros(size(x_cm2));
Af=zeros(sum(nr==1),1);
Isum=zeros(sum(nr==1),1);
Ir=zeros(sum(nr==1),Nrp);
g=zeros(sum(nr==1),Nrp);
i=0;

t=t(nr==1)';    % updates the time vector

for n_t=1:length(nr)
    % n_t==0 if a frame is omitted from the analysis
    if nr(n_t)==1
        i=i+1;

        I=double(imread(file_in,istart+n_t-1))-Idark;   % Reads the image
        I=reshape(I(BWframe==1),[length(y),length(x)]);
        I2=I./Ipre; % The relative intensity

        if i==1
            x_cm2(i)=x_cm;
            y_cm2(i)=y_cm;
            
            I20=I2;
            
            % Determines I2(inf) using values at the edge of the image
            nx=round(length(x)*0.1);
            ny=round(length(y)*0.1);
            I2s=sum(sum(I2(1:ny,:)))+sum(sum(I2(end:end-ny,:)))+...
                sum(sum(I2(:,1:nx)))+sum(sum(I2(:,end:end-nx)));
            I2n=sum(sum(~isnan(I2(1:ny,:))))+sum(sum(~isnan(I2(end:end-ny,:))))+...
                sum(sum(~isnan(I2(:,1:nx))))+sum(sum(~isnan(I2(:,end:end-nx))));
            I2_inf0=I2s/I2n;
        else
            x_cm2(i)=x_cm2(i-1);
            y_cm2(i)=y_cm2(i-1);
        end
      
        if i>2
            x_cm2(i)=(x_cm2(i-1)-x_cm2(i-2))/(t(i-1)-t(i-2))*(t(i)-t(i-1))+x_cm2(i-1);
            y_cm2(i)=(y_cm2(i-1)-y_cm2(i-2))/(t(i-1)-t(i-2))*(t(i)-t(i-1))+y_cm2(i-1);
        end
               
        % Take away the immobile fraction
        if trc=='y'
            gamma_g=0;   % immobile fraction
            I2=I2+(I2_inf0-I20)*gamma_g;
        end
         

        % Determines I2(inf) using values at the edge of the image
        nx=round(length(x)*0.1);
        ny=round(length(y)*0.1);
        I2s=sum(sum(I2(1:ny,:)))+sum(sum(I2(end:end-ny,:)))+...
            sum(sum(I2(:,1:nx)))+sum(sum(I2(:,end:end-nx)));
        I2n=sum(sum(~isnan(I2(1:ny,:))))+sum(sum(~isnan(I2(end:end-ny,:))))+...
            sum(sum(~isnan(I2(:,1:nx))))+sum(sum(~isnan(I2(:,end:end-nx))));
        I2_inf=I2s/I2n;
        
        % If trc=='y' then the program tracks (and moves) the center of mass
        % for each frame
        Iin=(I2_inf-I2);
        if trc=='y'
            %for j=1:3
                Rr2=(Xx-x_cm2(i)).^2+(Yy-y_cm2(i)).^2;
                BW=Rr2<(wid^2);
                cm_w=Iin(BW==1);
                Xx_w=Xx(BW==1);
                Yy_w=Yy(BW==1);

                Xv=Xx(BW);
                Yv=Yy(BW);
                Iv=I2(BW)/mean(I2(:));
                
                % Calculates the centre of mass (x_cm,y_cm)
                pcm=lsqnonlin(@(p) gauss_fkn(p,Xv,Yv,Iv),pcm0,...
                    [0,-Inf,min(Xv),min(Yv),0],[Inf,Inf,max(Xv),max(Yv),Inf],options);
                
                x_cm2(i)=pcm(3);
                y_cm2(i)=pcm(4);
                pcm0=pcm;
                
                %x_cm2(i)=sum(Xx_w.*cm_w)/sum(cm_w);
                %y_cm2(i)=sum(Yy_w.*cm_w)/sum(cm_w);
            %end   
        else      
            Rr2=(Xx-x_cm2(i)).^2+(Yy-y_cm2(i)).^2;
            BW=Rr2<(wid^2);
            cm_w=Iin(BW==1);
            Xx_w=Xx(BW==1);
            Yy_w=Yy(BW==1);
            x_cm2(i)=x_cm;
            y_cm2(i)=y_cm;           
        end
        
        % Calculates the moments (Ixx and Iyy) and product (Ixy) of inertia for the bleached spot
        Ixx(i)=sum((Yy_w-y_cm2(i)).^2.*cm_w);
        Iyy(i)=sum((Xx_w-x_cm2(i)).^2.*cm_w);
        Ixy(i)=sum((Xx_w-x_cm2(i)).*(Yy_w-y_cm2(i)).*cm_w);

        % Performs the radial averaging of I2 around the centre of mass
        Ir(i,:)=I_radial(Nrp,x_cm2(i),y_cm2(i),Xx,Yy,I2,Rp_max);

        % Determines Ir(R,0) using values where r>0.9R
        if i==1
            Ir0=mean(Ir(1,r>(max(r)*0.9)));
        end

        % Determines the integral of Ir at r<R
        Isum(i)=2*pi*Ir(i,1)*r(1)*r(1);
        for j=2:length(r)
            Isum(i)=Isum(i)+2*pi*(Ir(i,j)*r(j)+Ir(i,j-1)*r(j-1))/2*(r(j)-r(j-1));
        end

        % Determines starting parameters used to fit the tail of Ir
        yf=Ir(i,:)/Ir0;
        xf=r;
        pos=1;
        yf0=yf(pos);
        while abs(1-yf0)>abs(1-yf(1))*exp(-1) && pos<length(yf)
            pos=pos+1;
            yf0=yf(pos);
        end
        w_g0=sqrt(xf(pos)^2-xf(1)^2);

        yf=Ir(i,r>max(r)*1/2)/Ir0;
        xf=r(r>max(r)*1/2);
        pos=1;
        yf0=yf(pos);
        while abs(1-yf0)>abs(1-yf(1))*exp(-1) && pos<length(yf)
            pos=pos+1;
            yf0=yf(pos);
        end
        w_g=max(sqrt(xf(pos)^2-xf(1)^2),w_g0);

        % Makes a curve fit to the tail of Ir (for r>R/2) for each frame
        p0=[sign(1-yf(1))*min(abs(1-Ir(i,1)/Ir0),abs(1-yf(1))*exp(xf(1)^2/w_g^2)),w_g];
        p0_2=p0;
        if i>1
            p0_2=p;
        end
        
        [p1,sse1]=lsqcurvefit(@(p,r) fkn_gauss(p,r,Isum(i),Isum(1),Ir0),p0,xf,yf*Ir0,...
            [-Inf,0],[Inf,Inf],options);

        [p2,sse2]=lsqcurvefit(@(p,r) fkn_gauss(p,r,Isum(i),Isum(1),Ir0),p0_2,xf,yf*Ir0,...
            [-Inf,0],[Inf,Inf],options);

        if sse1<sse2
            p=p1;
        else
            p=p2;
        end
        
        wf(i)=p(2); % w(t)
        if i==1
            p(1)=0;
        end
        Af(i)=p(1); % A(t)
        [yf,beta]=fkn_gauss(p,r,Isum(i),Isum(1),Ir0);    % y = Ir|r>R

        % Compensates for temporal variations using beta and for
        % a net influx of molecules with y
        g(i,:)=(yf-Ir(i,:))/beta;
               
        if mod(i,10)==0 || i==1 || i==sum(nr==1)
            figure(2)
            ax1 = gca;
            set(ax1,'FontSize',12)
            plot(r,Ir(i,:)/beta,'b.',r(r>(max(r)*1/2)),yf(r>(max(r)*1/2))/beta,'r-','LineWidth',2)
            xlabel('r [\mum]','FontSize',12)
            ylabel('Intensity [a.u.]','FontSize',12)
            title(['Frame ',num2str(i),' out of ',num2str(sum(nr==1))],'FontSize',12)
            xlim([min(r),max(r)])
            pause(0.1)
        end

    end
end

% Calculates the Hankel transform and performs the curve fits yielding D1,
% D2, gamma2 and gamma0
[k,D1k,D1,D2,gamma2,gamma0,v_conv,a_conv]=...
    hankel_diff(r,t,R,g,wf,Af,Nrp,trc,x_cm2,y_cm2,dx,Ixx,Iyy,Ixy);

end

function dIv=gauss_fkn(pcm,Xv,Yv,Iv)
    dIv=pcm(1)*(1-pcm(2)*exp(-((Xv-pcm(3)).^2+(Yv-pcm(4)).^2)/pcm(5)^2))-Iv;
end
