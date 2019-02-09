///////////////////////////////////////////////////
// PSM - Point Shift Measurement
// Scilab programm
// For BOS technique
// requires installation IPCV atom
// first version - 02-2019 
// Jury Barinov


///////////////////////////////////////////////////
// Image_1 preprocessing
// remove background, increase contrast
function FrO=Im_pr1(Fr)
    disp("Preprocessing File #1");
    //to marix
    FrO = im2double(Fr)
    // negativity
    FrO=imcomplement(FrO);
    [y,x,z]=size(FrO)
    if z==3 then
        disp("File #1 is color")
        // grayscale
        FrO = rgb2gray(FrO)
    end    
    // Setting is bright. and contrast
    FrO=imadjust(FrO,[0.2 1],[0 1]); // The options depend on the image quality
    // Noise filtering
    filter = fspecial('average',7)  // 
    FrO = imfilter(FrO, filter);
    // Processing by pieces
    // Similar to adaptive level
    FrO=AdThreshold(FrO,3,5)
endfunction

///////////////////////////////////////////////////
//Image_2 preprocessing
// remove background, increase contrast
function FrO=Im_pr2(Fr)
    disp("Preprocessing File #2");
    //to marix
    FrO = im2double(Fr)
    // negativity
    FrO=imcomplement(FrO);
    [y,x,z]=size(FrO)
    if z==3 then
        disp("File #2 is color")
        // grayscale
        FrO = rgb2gray(FrO)
    end
    // Setting is bright. and contrast
    FrO=imadjust(FrO,[0.05 1],[0 1]); // The options depend on the image quality
    // Noise filtering
    filter = fspecial('average',7)  // хорошо
    FrO = imfilter(FrO, filter);    
    // Processing by pieces
    // Similar to adaptive level
    FrO=AdThreshold(FrO,7,4)    
endfunction

///////////////
// Binarization by parts
// I-input im.
// k - number of pieces of partition
function [B]=AdThreshold(I,ky,kx)
    //
    [yp,xp]=size(I)
    dy = fix(yp/ky)
    dx = fix(xp/kx)
    dyy = yp - ky*dy
    dxx = xp - kx*dx
    B=zeros(I)
    for j=1:ky
        for i=1:kx
            y=(j-1)*dy+1:j*dy
            x=(i-1)*dx+1:i*dx
            th = imgraythresh( I( y,x ) );
            B(y,x) = im2bw(I(y,x),th);
        end 
    end     
    for i=1:kx
        y=ky*dy+1:ky*dy+dyy
        x=(i-1)*dx+1:i*dx
        th = imgraythresh( I( y,x ) );
        B(y,x) = im2bw(I(y,x),th);

    end      
    for j=1:ky-1
        y=(j-1)*dyy+1:j*dyy
        x=kx*dx+1:kx*dx+dxx
        th = imgraythresh( I( y,x ) );
        B(y,x) = im2bw(I(y,x),th);
    end
endfunction

///////////////////////////////////////////////////
// find points 
// Fr1 - background screen
// Fr2 - distorted image
function [ctr1, ctr2, BB1,BB2,A1,A2]=Detect_point(Fr1, Fr2)
    disp("Detect blobs...")
    // edge search
    disp("For #1")
    S_edge1 = edge(Fr1, 'fftderiv',0.3); // the ratio depends on the image quality
    [A1,BB1,ctr1]=Fr_Bkg(S_edge1)
    // edge search
    disp("For #2")
    S_edge2 = edge(Fr2, 'fftderiv',0.05); // the ratio depends on the image quality
    [A2,BB2,ctr2]=Fr_Mv(S_edge2)
    disp("Time="+string(toc())+"s")
endfunction

//////////////////////////
// point search function on the background screen
// Fr-background screen photo
// BB-selected areas
// cntr0-output array of point centers in the selected area
function [A,BB,cntr0]=Fr_Bkg(Fr)
    disp("Time="+string(toc())+"s")
    [S_labeled1,n1] = imlabel(Fr);
    disp("Found blobs="+string(n1))
    disp("Time="+string(toc())+"s")
    // find blobs
    [A, BB,cntr0] = _imblobprop(S_labeled1)

endfunction

///////////////////////////////////////////////////
// search for offset points on distorted image
// Fr - photo of the background screen with offset points; 
// BB-selected areas
// cntr0-output array of point centers in the selected area
function [A,BB,cntr0]=Fr_Mv(Fr)
    disp("Time="+string(toc())+"s")
    [S_labeled2,n2] = imlabel(Fr);
    disp("Found blobs="+string(n2))
    disp("Time="+string(toc())+"s")
    // find blobs
    [A, BB,cntr0] = _imblobprop(S_labeled2)
endfunction

///////////////////////////////////////////////////
// find neighbor points
// searches for points close to the background screen points
// cntr1-background screen points
// cntr2-shift points 
// rd is the search radius
function[su,dc,Nn]= find_poins(cntr1, cntr2,A1,A2,rd)
    disp("F7")
    dc=[]
    h = length(cntr1(1,:))
    Nn = zeros(h,10)
    N=0
    for i=1:h
        n=find( ( (cntr1(1,i)-cntr2(1,:))^2 + (cntr1(2,i)-cntr2(2,:))^2 )<=rd^2 ) // search within radius
        if ~isempty(n)then 
            dc(2,i)=mean(cntr1(2,i)-cntr2(2,n))// the average distance in x
            dc(1,i)=mean(cntr1(1,i)-cntr2(1,n))// the average distance in y
            su(1,i)=-(dc(1,i)-cntr1(1,i))  // x shift
            su(2,i)=-(dc(2,i)-cntr1(2,i))  // y shift          
            Nn(i,1:length(n)) = n
            N=N+length(n)
        else 
            su(1,i)=cntr1(1,i)
            su(2,i)=cntr1(2,i)
        end
    end
    if N < 0.1*h then // if images are very different
        disp("Files are very different")
    else
        disp("Blobs for calc.="+string(N))
    end
endfunction    

///////////
// Removing the wrong points
function [A,cntr,BB]=sort(A,cntr,BB)
    disp("Sort")
    mx=max(A)
    mn=min(A)
    ma=mean(A)
    if mx>=2*ma then
        n=find(A>mx*0.7)
        A(n)=[]
        cntr(:,n)=[]
        BB(:,n)=[]
        disp("dell(1)="+string(length(n)))
    end
    n=find(A<0.6*mean(A))
    disp("dell(2)="+string(length(n)))
    cntr(:,n)=[]
    BB(:,n)=[]
endfunction

///////////
// Other IPCV function "imblobprop" implementation
// Works about three times faster
function [A, BB, ctr] = _imblobprop(imin)
    n = double(max(imin));
    A = zeros(1,n);
    BB = zeros(4,n);
    ctr = zeros(2,n);    
    for cnt = 1:n
        [r,c] = find(imin==cnt);
        A(cnt) = length(r)
        minx = min(c);
        miny = min(r);
        maxx = max(c);
        maxy = max(r);
        w = maxx-minx; // +1
        h = maxy-miny; // +1
        BB(:,cnt) = [minx; miny; w; h];
        avgx = mean(c);
        avgy = mean(r);
        ctr(:,cnt) = [avgx; avgy];
    end
endfunction

///////////////
// Finding the shift
function [i,xm,ym]=shift(BB1,BB2,Nn)
    i=find(Nn(:,1)>0)    
    xm1=BB1(1,i)-BB2(1,Nn(i,1))
    xm=2*xm1+(BB1(3,i)-BB2(3,Nn(i,1))) 
    ym1=BB1(2,i)-BB2(2,Nn(i,1))
    ym=2*ym1+(BB1(4,i)-BB2(4,Nn(i,1))) 
endfunction

///////////
// Averaging
function [out,X,Y]=avr(cntr,mv)
    //
    dx=5  // averaging range
    //  
    dy=round(max(cntr(2,:))/3) // // averaging range to y=3
    X=[round(min(cntr(1,:))):dx:round(max(cntr(1,:)))]
    Y=[round(min(cntr(2,:))):dy:round(max(cntr(2,:)))]    
    out=double(zeros(length(X),length(Y)))
    for x=1:length(X)
        nx=find((cntr(1,:)<(X(x)+0.5*dx)) & (cntr(1,:)>(X(x)-0.5*dx)))
        //nx=i(ni)
        //disp(nx)
        if ~isempty(nx)then 
            for y=1:length(Y)
                ny=find((cntr(2,nx)<(Y(y)+0.5*dy)) & (cntr(2,nx)>(Y(y)-0.5*dy)))
                //disp(ny)
                if ~isempty(ny)then
                    out(x,y)=mean(mv(nx(ny)))
                else
                    out(x,y)=%nan
                end
            end
        else 
            out(x,:)=%nan
        end
    end  
endfunction


////////////////////////////
//        Main
////////////////////////////

// To insert a file path
filename_1="F:/TMP/file1.png"//"Image file path" // background image
filename_2="F:/TMP/file2.png"//"Image file path" // distorted image

S1=imread(filename_1)
S2=imread(filename_2)
        
rd=15 // search radius


S12=Im_pr1(S1)
S22=Im_pr2(S2)

tic();

[ys,xs]=size(S1)
f1 = scf();
f1.figure_size= [xs,ys];
subplot(1,2,1)   
imshow(S12);
xtitle("File #1")
//
subplot(1,2,2)
imshow(S22);
xtitle("File #2")
//

disp("Time="+string(toc())+"s")
//
[ctr0, ctr1, BB1, BB2,A1,A2]=Detect_point(S12, S22)
disp("Time="+string(toc())+"s")
//
[A1,ctr0,BB1]=sort(A1,ctr0,BB1)
[A2,ctr1,BB2]=sort(A2,ctr1,BB2)
//
[su,dc,Nn]= find_poins(ctr0, ctr1, A1,A2,rd)
disp("Time="+string(toc())+"s")
//
[i,xm,ym]=shift(BB1,BB2,Nn)

f2 = scf();
gcf().color_map =  hotcolormap(64); 
subplot(1,2,1)
scatter3(ctr0(1,i),ys-ctr0(2,i)+1,(xm),5,xm)
xtitle("Shift X")
subplot(1,2,2)
scatter3(ctr0(1,i),ys-ctr0(2,i)+1,(ym),5,ym)
xtitle("Shift Y")
//
[out,X,Y]=avr(ctr0(:,i),xm)
//

f3 = scf();
plot(X,out(:,:),"b." )
xgrid(color("grey"));
xtitle("All shift points X")
//

// If necessary, delete comments
///////////////
// save to file
//yy=out(:,1)' // 1 or 2 ...
//Dell Nan
//k = find(~isnan(out(:,1))) // 1 or 2 ...
//yy = out(k,1) // 1 or 2 ...
//Mx=[X(k)',yy]
// saves only the x shift.

//savepath="path_to_save"
//filesave="name_file"
//csvWrite(Mx, savepath+filesave+".txt",ascii(9));
//    

