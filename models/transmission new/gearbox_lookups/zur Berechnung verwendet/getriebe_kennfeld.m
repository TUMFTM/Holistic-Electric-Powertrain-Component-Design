% File for computing a lookup-table with Gao Gearbox Model

tic;
loss_matrix = zeros(42,21,2);

for t=1:2 %temperature 
    temp = 25*t;
    for i = 0:20 %rpm
        rpm_ind = 75*i+10
        for j = 0:41 %torque
            if j<21
                torque_ind = -(2200-110*j+20)
            else
            torque_ind = (110*(j-21)+20)
            end
    
            system(strcat('python Transmissionpy_Kennfeld.py --tp'," ",sprintf('%.0f',3)," ",'--temp'," ",sprintf('%f',temp)," ",'--torque'," ",sprintf('%f',torque_ind)," ",'--rpm '," ",sprintf('%f',rpm_ind)));
            
            fid = fopen('transmission_output.txt','r');
            tline = fgets(fid);    
            data = sscanf(tline,'%f');
            powerloss = data;
            fclose(fid);
    
            loss_matrix(j+1,i+1,t) = powerloss;
    
        end
    end
end
toc