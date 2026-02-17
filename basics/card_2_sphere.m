function [ A_th , A_ph , A_R ] = card_2_sphere( A , TH , PH )
    a_x = A( : , : , 1 ) ;
    a_y = A( : , : , 2 ) ;
    a_z = A( : , : , 3 ) ;
    
    A_th = a_x.*cos(TH).*cos(PH) + a_y.*cos(TH).*sin(PH) - a_z.*sin(TH) ;
    A_ph = - a_x.*sin(PH) + a_y.*cos(PH) ;
    A_R = a_x.*sin(TH).*cos(PH) + a_y.*sin(TH).*sin(PH) + a_z.*cos(TH) ;

end

