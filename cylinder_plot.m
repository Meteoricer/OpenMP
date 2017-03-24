%temp_data=load('output - Copy.txt');
fid=fopen('output.txt');
col_num=600;
row_num=800;
fig=figure(1);
set(fig,'Position',[100,100,900,900]);
filename='cylinder_plot_with_ring_formed.gif';
data=zeros(col_num,col_num);
line_number=0;
R=600/2/pi*5;
mark=0;
%temp=zeros(col_num+2:row_num+2:3);
while 1
nextline = fgetl(fid); %read a line 
    if ~ischar(nextline)
        fclose(fid);
        break;
    else
        line_number = line_number + 1
        if(line_number>1300)
            temp=ones(col_num+2,row_num+2,3)*255;
            count=0;
            temp_data=str2num(nextline);
            temp_x=zeros(length(temp_data)-1,1);
            temp_y=zeros(length(temp_data)-1,1);
            for i=2:2:length(temp_data)
            %templa(i)=temp_data(i);
                count=count+temp_data(i);
                col=rem(count,row_num)+1;
                row=floor(count/row_num)+1;
                %data(row,col)=temp_data(i);
                switch temp_data(i+1)
                    case 1
                         temp_y(i-1)=row;
                         temp_x(i-1)=col;
                         %orient_data(f,1)=orient_data(f,1)+1;


                    case 2
                         temp_y(i-1)=row;
                         temp_x(i-1)=col;
                         %orient_data(f,2)=orient_data(f,2)+1;

                    case 3
                         temp_y(i-1)=row;
                         temp_x(i-1)=col;
                    case 4
                         temp_y(i-1)=row;
                         temp_x(i-1)=col;


                    case 5
                         temp_y(i-1)=row;
                         temp_x(i-1)=col;
                      
                    
                end    
            
            end
        
            theta=2*pi*temp_y/600;
            x=R*cos(theta);
            y=R*sin(theta);
            z=temp_x*5;





            scatter3(z,x,y,1,'filled');
            
            axis([0 4000 -2000 2000 -2000 2000]);
            frame = getframe(fig);
            im = frame2im(frame);
            [A,map]=rgb2ind(im,256);
            if mark==0
                imwrite(A,map,filename,'gif','LoopCount',inf,'DelayTime',0.1);
                mark=1;
            else 
                imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',0.1);
            end
            
        end
        

       
    end
end 
test_data=zeros(3000,6);
% for f=1:1:3000
%     for i=2:2:length(temp_data(f))
%         temp_data(f,i)
%     end
% end

