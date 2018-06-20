function  [struct1,struct1_Namelist,struct2,struct2_Namelist] = common_str_final_AF(input1,input1_Namelist,input2,input2_Namelist)

%input1 = protein;
%input1_AtomName = protein_AtomName;
%input1_resName = protein_resName;
%input2 = Apo_protein;
%input2_AtomName = Apo_protein_AtomName;
%input2_resName = Apo_protein_resName;


input1(:,10)=0;
input2(:,10)=0;

i=1;
j=1;
struct1=[];
struct2=[];
struct1_Namelist = [];
struct2_Namelist = [];

while i<=size(input1,1)
    itt=0;
    
    while j<=size(input2,1)
        if i+50<=size(input1,1) && j+50<=size(input2,1)
            

            if all(input1(i:i+50,5)==input2(j:j+50,5)) == 1 && all(input1(i:i+50,10)==0) == 1 && all(input2(j:j+50,10)==0) == 1
                
                for k =1:size(input2,1)-j
                    if i+k <=size(input1,1) && j+k <=size(input2,1)
                        
                
                        if all(input1(i:i+k,5)==input2(j:j+k,5)) == 1 && (i+k < size(input1,1) && j+k < size(input2,1)) && all(input1(i:i+k,10)==0) == 1 && all(input2(j:j+k,10)==0) == 1
                            continue
                        elseif all(input1(i:i+k,5)==input2(j:j+k,5)) == 1 && (i+k < size(input1,1) && j+k < size(input2,1)) && (all(input1(i:i+k,10)==0) ~= 1 || all(input2(j:j+k,10)==0) ~= 1)
                            
                            struct1_Namelist(size(struct1,1)+1:size(struct1,1)+k,:)=input1_Namelist(i:i+k-1,:);
                            struct2_Namelist(size(struct2,1)+1:size(struct2,1)+k,:)=input2_Namelist(j:j+k-1,:);
                            
                            struct1(size(struct1,1)+1:size(struct1,1)+k,:)=input1(i:i+k-1,:);
                            struct2(size(struct2,1)+1:size(struct2,1)+k,:)=input2(j:j+k-1,:);
                            
                            input1(i:i+k-1,10)=1;
                            input2(j:j+k-1,10)=1;
                            i=i+k;
                            j=j+k;
                            itt=1;
                            break
                        elseif all(input1(i:i+k,5)==input2(j:j+k,5)) == 1 & (i+k == size(input1,1) | j+k == size(input2,1)) && all(input1(i:i+k,10)==0) == 1 && all(input2(j:j+k,10)==0) == 1
                            
                            struct1_Namelist(size(struct1,1)+1:size(struct1,1)+k+1,:)=input1_Namelist(i:i+k,:);
                            struct2_Namelist(size(struct2,1)+1:size(struct2,1)+k+1,:)=input2_Namelist(j:j+k,:);
                            struct1(size(struct1,1)+1:size(struct1,1)+k+1,:)=input1(i:i+k,:);
                            struct2(size(struct2,1)+1:size(struct2,1)+k+1,:)=input2(j:j+k,:);
                            input1(i:i+k,10)=1;
                            input2(j:j+k,10)=1;
                            i=i+k;
                            j=j+k;
                            itt=1;
                            break    
                        elseif all(input1(i:i+k,5)==input2(j:j+k,5)) == 1 & (i+k == size(input1,1) | j+k == size(input2,1)) && (all(input1(i:i+k,10)==0) ~= 1 || all(input2(j:j+k,10)==0) ~= 1)
                            
                            struct1_Namelist(size(struct1,1)+1:size(struct1,1)+k,:)=input1_Namelist(i:i+k-1,:);
                            struct2_Namelist(size(struct2,1)+1:size(struct2,1)+k,:)=input2_Namelist(j:j+k-1,:);
                            
                            struct1(size(struct1,1)+1:size(struct1,1)+k,:)=input1(i:i+k-1,:);
                            struct2(size(struct2,1)+1:size(struct2,1)+k,:)=input2(j:j+k-1,:);
                            input1(i:i+k-1,10)=1;
                            input2(j:j+k-1,10)=1;
                            i=i+k;
                            j=j+k;
                            itt=1;
                            break

                        elseif all(input1(i:i+k,5)==input2(j:j+k,5)) ==0 && (i+k < size(input1,1) && j+k < size(input2,1)) && all(input1(i:i+k-1,10)==0) == 1 && all(input2(j:j+k-1,10)==0) == 1
                            
                            struct1_Namelist(size(struct1,1)+1:size(struct1,1)+k,:)=input1_Namelist(i:i+k-1,:);
                            struct2_Namelist(size(struct2,1)+1:size(struct2,1)+k,:)=input2_Namelist(j:j+k-1,:);
                            
                            struct1(size(struct1,1)+1:size(struct1,1)+k,:)=input1(i:i+k-1,:);
                            struct2(size(struct2,1)+1:size(struct2,1)+k,:)=input2(j:j+k-1,:);
                            
                            input1(i:i+k-1,10)=1;
                            input2(j:j+k-1,10)=1;
                            i=i+k;
                            j=j+k;
                            itt=1;
                            break
                        elseif all(input1(i:i+k,5)==input2(j:j+k,5)) ==0 && (i+k < size(input1,1) && j+k < size(input2,1)) && all(input1(i:i+k-1,10)==0) ~= 1 && all(input2(j:j+k-1,10)==0) ~= 1
                            
                            
                            struct1_Namelist(size(struct1,1)+1:size(struct1,1)+k,:)=input1_Namelist(i:i+k-1,:);
                            struct2_Namelist(size(struct2,1)+1:size(struct2,1)+k,:)=input2_Namelist(j:j+k-1,:);
                            
                            
                            struct1(size(struct1,1)+1:size(struct1,1)+k,:)=input1(i:i+k-1,:);
                            struct2(size(struct2,1)+1:size(struct2,1)+k,:)=input2(j:j+k-1,:);
                            input1(i:i+k-1,10)=1;
                            input2(j:j+k-1,10)=1;
                            i=i+k;
                            j=j+k;
                            itt=1;
                            break
                            
                            
                        elseif all(input1(i:i+k,5)==input2(j:j+k,5)) ==0 && (i+k == size(input1,1) || j+k == size(input2,1)) && all(input1(i:i+k-1,10)==0) == 1 && all(input2(j:j+k-1,10)==0) == 1
                            
                            struct1_Namelist(size(struct1,1)+1:size(struct1,1)+k,:)=input1_Namelist(i:i+k-1,:);
                            struct2_Namelist(size(struct2,1)+1:size(struct2,1)+k,:)=input2_Namelist(j:j+k-1,:);
                            
                            struct1(size(struct1,1)+1:size(struct1,1)+k,:)=input1(i:i+k-1,:);
                            struct2(size(struct2,1)+1:size(struct2,1)+k,:)=input2(j:j+k-1,:);
                            input1(i:i+k-1,10)=1;
                            input2(j:j+k-1,10)=1;
                            i=i+k;
                            j=j+k;
                            itt=1;
                            break
                        elseif all(input1(i:i+k,5)==input2(j:j+k,5)) ==0 && (i+k == size(input1,1) || j+k == size(input2,1)) && all(input1(i:i+k-1,10)==0) ~= 1 && all(input2(j:j+k-1,10)==0) ~= 1
                            
                            struct1_Namelist(size(struct1,1)+1:size(struct1,1)+k,:)=input1_Namelist(i:i+k-1,:);
                            struct2_Namelist(size(struct2,1)+1:size(struct2,1)+k,:)=input2_Namelist(j:j+k-1,:);
                            
                            struct1(size(struct1,1)+1:size(struct1,1)+k,:)=input1(i:i+k-1,:);
                            struct2(size(struct2,1)+1:size(struct2,1)+k,:)=input2(j:j+k-1,:);
                            input1(i:i+k-1,10)=1;
                            input2(j:j+k-1,10)=1;
                            i=i+k;
                            j=j+k;
                            itt=1;
                            break   
                            
                        
                        end
                    elseif i+k > size(input1,1) || j+k > size(input2,1)
                        j=j+1;
                        break
                    end
                end
                    
            elseif all(input1(i:i+50,5)==input2(j:j+50,5)) ~= 1 || all(input1(i:i+50,10)==0) ~= 1 || all(input2(j:j+50,10)==0) ~= 1
                j=j+1;
            end
        elseif i+50>size(input1,1) || j+50>size(input2,1)
            j = size(input2,1)+1;
            
        end
    end
    if itt == 0
        i=i+1;
        j=1;
    end
end
struct1(:,10) = 1:size(struct1,1);  
struct2(:,10) = 1:size(struct2,1);                         
                        
                    
