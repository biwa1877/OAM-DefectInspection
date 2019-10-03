function U = LaguerreGaussianBeam_v00(p, l, w0, lambda, z, r_edge)
    
    if isempty(p)
        p = 0;                  % Degree of LG mode
    end
    
    if isempty(l)
        l = 1;                  % Order of LG mode
    end
    
    if isempty(w0)
        w0 = 15;               % Beam waist
    end
    
    if isempty(lambda)
        lambda = 532.0e-9;
    end
    k = 2*pi / lambda;      % Wavenumber of light

    zR = k*w0^2/2;      % Calculate the Rayleigh range

    % Setup the cartesian grid for the plot at plane z
    if isempty(z)
        z = 0.0;
    end
    
    % [xx, yy] = meshgrid(linspace(-8, 7), linspace(-8, 7));
    [xx, yy] = meshgrid(-128:128);

    % make a hard edge
    if isempty(r_edge)
        r_edge = 50;
    end
    
    hardEdge = ( (xx.^2 + yy.^2) <= r_edge^2 );

    % Calculate the cylindrical coordinates
    [phi, r] = cart2pol(xx, yy);

    U00 = 1/(1 + 1i*z/zR) .* exp(-r.^2/w0^2./(1 + 1i*z/zR));
    w = w0 * sqrt(1 + z.^2/zR^2);
    R = sqrt(2)*r./w;

    % Lpl from OT toolbox (Nieminen et al., 2004)
    Lpl = nchoosek(p+l,p) * ones(size(R));   % x = R(r, z).^2
    for m = 1:p
        Lpl = Lpl + (-1)^m/factorial(m) * nchoosek(p+l,p-m) * R.^(2*m);
    end

    U = U00.*R.^l.*Lpl.*exp(1i*l*phi).*exp(-1i*(2*p + l + 1)*atan(z/zR));

    % give the beam a hard edge
    U = U .* hardEdge;

%     figure(101);
%     subplot(1, 2, 1);
%     % imagesc(log(abs(U).^2));axis image;
%     imagesc((abs(U).^2));axis image;
%     title('Intensity');
%     subplot(1, 2, 2);
%     imagesc(angle(U));axis image;
%     title('Phase');

end