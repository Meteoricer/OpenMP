%temp_data=load('output - Copy.txt');
fid=fopen('output.txt');
col_num=600;
row_num=800;
figure(1);
filename='test.gif';
data=zeros(col_num,col_num);
line_number=0;
r1=300;     %????
r2=r1+50;    %????

%temp=zeros(col_num+2:row_num+2:3);
while 1
nextline = fgetl(fid); %read a line 
    if ~ischar(nextline)
        fclose(fid);
        break;
    else
        line_number = line_number + 1
        temp=ones(col_num,row_num,3)*255;
        img=ones(row_num,col_num)*255;
        image=ones(2*r2,2*r2,3)*255;
        count=0;
        temp_data=str2num(nextline);
        for i=2:2:length(temp_data)
            %templa(i)=temp_data(i);
            count=count+temp_data(i);
            col=rem(count,row_num)+1;
            row=floor(count/row_num)+1;
            %data(row,col)=temp_data(i);
            if (col>left && col<right)
                 switch temp_data(i+1)
                
                
                
                
                    case 1
                         temp(row,col-left,1)=255;
                         temp(row,col-left,2)=0;
                         temp(row,col-left,3)=0;
                         %orient_data(f,1)=orient_data(f,1)+1;


                    case 2
                         temp(row,col-left,1)=0;
                         temp(row,col-left,2)=0;
                         temp(row,col-left,3)=255;
                         %orient_data(f,2)=orient_data(f,2)+1;

                    case 3
                         temp(row,col-left,1)=0;
                         temp(row,col-left,2)=0;
                         temp(row,col-left,3)=0;
                    case 4
                         temp(row,col-left,1)=255;
                         temp(row,col-left,2)=0;
                         temp(row,col-left,3)=0;


                    case 5
                         temp(row,col-left,1)=0;
                         temp(row,col-left,2)=0;
                         temp(row,col-left,3)=255;


                end            
            end
               
            
        end
        
        
        for k=1:1:3
            img=transpose(temp(:,:,k));
            %imshow(img);
            [m,n]=size(img);

            
            imgn=zeros(2*r2,2*r2);
            [re_m,re_n]=size(imgn);
            for y=1:re_m
                for x=1:re_n
                    dis_x=x-re_n/2;
                    dis_y=y-re_m/2;

                    l=sqrt(dis_x^2+dis_y^2);
                    if l<=r2 && l>=r1
                        theta=0;
                        if y>re_m/2
                            theta=atan2(dis_y,dis_x);
                        end
                        if y<re_m/2
                            theta=pi+atan2(-dis_y,-dis_x);
                        end            
                        if y==re_m/2
                            theta=atan2(dis_y,dis_x)+0.0001;
                        end

                        xx=ceil(n*theta/(2*pi));
                        yy=ceil(l-r1);
                        if yy>=1 && yy<=m && xx>=1 && xx<=n
                            imgn(y,x)=img(yy,xx);
                        end
                    end
                end
            end
            image(:,:,k)=imgn;
            
        end
        

        
        
        
        
        
        imshow(image);
        frame = getframe;
        im = frame2im(frame);
        [A,map]=rgb2ind(im,256);
        if line_number==1
            imwrite(A,map,filename,'gif','LoopCount',inf,'DelayTime',0.1);
        else 
            imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',0.1);
        end

       
    end
end 
test_data=zeros(3000,6);
% for f=1:1:3000
%     for i=2:2:length(temp_data(f))
%         temp_data(f,i)
%     end
% end

