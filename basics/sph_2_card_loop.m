function [ A ] = sph_2_card_loop( a_th , a_ph , a_R , TH , PH )

    A(:,:,1) = a_th.*cos(TH).*cos(PH) - a_ph.*sin(PH) + a_R.*sin(TH).*cos(PH) ;
    A(:,:,2) = a_th.*cos(TH).*sin(PH) + a_ph.*cos(PH) + a_R.*sin(TH).*sin(PH) ;
    A(:,:,3) = - a_th.*sin(TH) + a_R.*cos(TH) ;

end

