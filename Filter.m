new_test_data=test_data(begin:4324,:);
max_num=max(test_data);
min_num=mean(test_data(begin:4324,:));
[width,length_data]=size(length_of_ring);
for i=1:num_boo
    new_test_data(:,i)=new_test_data(:,i)-0;
end
critical=0;



for i=1:2325
     
         
         for k=1:1:num_boo
             if(new_test_data(i,k)<0)
                 new_test_data(i,k)=0;
             else
                 %new_test_data(i,k)=1;
             end
         end
         
     

end
