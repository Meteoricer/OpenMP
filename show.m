data=load('output - Copy.txt');
col_num=600;
row_num=800;
figure(1);
filename='test.gif';
temp=zeros(col_num+2:row_num+2:3);
orient_data=zeros(1800,2);
test_data=zeros(1800);
%temp=zeros(col_num,row_num,1,600);
%M=moviein(20);
i=1;
for j=1:1:row_num+2
    temp(i,j,1)=0;
    temp(i,j,2)=0;
    temp(i,j,3)=0;
end
i=col_num+2;
for j=1:1:row_num+2
    temp(i,j,1)=0;
    temp(i,j,2)=0;
    temp(i,j,3)=0;
end
for f=1:1:1800
    %temp=zeros(col_num+2,row_num+2);
   
    for i=2:1:col_num+1
        for j=1:1:row_num+1
%            if (j==1 || j==row_num+2)
%                temp(i,j)=0;
%            else
%                if data(f,(i-2)*row_num+j)==1%T
%                    temp(i,j)
%                temp(i,j)=~data(f,(i-2)*row_num+j);
%            end
            switch data(f,(i-2)*row_num+j)
                case 0
                     temp(i,j,1)=255;
                     temp(i,j,2)=255;
                     temp(i,j,3)=255;
                case 1
                     temp(i,j,1)=255;
                     temp(i,j,2)=0;
                     temp(i,j,3)=0;
                     orient_data(f,1)=orient_data(f,1)+1;
                     if (j>130 && j<280)&&(i>200&&i<400)
                         test_data(f)=test_data(f)+1;
                     end
                case 2
                     temp(i,j,1)=0;
                     temp(i,j,2)=0;
                     temp(i,j,3)=255;
                     orient_data(f,2)=orient_data(f,2)+1;
                     if (j>130 && j<280)&&(i>200&&i<400)
                         test_data(f)=test_data(f)+1;
                     end
                case 3
                     temp(i,j,1)=0;
                     temp(i,j,2)=0;
                     temp(i,j,3)=0;
                    
            end    
            
        end

    end
    imshow(temp);
    frame = getframe;
    im = frame2im(frame);
    [A,map]=rgb2ind(im,256);
    if f==1
        imwrite(A,map,filename,'gif','LoopCount',inf,'DelayTime',0.1);
    else 
        imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',0.1);
    end
   
    %M(:,f)=getframe;
end

%imwrite(temp,filename,'gif')