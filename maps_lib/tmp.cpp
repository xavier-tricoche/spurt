double rotation_angle_linear_predictor(const nvis::vec2& x0, const nvis::vec2& x1,
                                       // v0 and v1 are the values associated with x0 and x1...
                                       const nvis::vec2& v0, const nvis::vec2& v1,
                                       // ...provided their respective flag valid? is true
                                       bool valid0, bool valid1,
                                       const integrator_type& map, unsigned int period,
                                       double dtheta, double dx,
                                       const metric_type& metric,
                                       unsigned int depth,
                                       std::list<step_type>& steps,
                                       bool conservative = false) const
{
    // needed for internal consistency
    const double large_angle = 0.75 * M_PI;

    nvis::vec2 __v0, __v1, dummy;
    __v0 = (valid0 ? v0 : evaluate_step(x0, map, period, metric).first);
    __v1 = (valid1 ? v1 : evaluate_step(x1, map, period, metric).first);

    // easy cases first
    double theta = angle(__v0, __v1);
    if (fabs(theta) < std::max(large_angle, dtheta)) return theta;
    else if (nvis::norm(x1 - x0) < dx) {
        if (conservative) {
            throw degenerate_point_exception(0.5*(x0 + x1));
        }
        return theta;
    }

    // use linear model to pick next sample
    // solve for v(x).v0 = 0
    double v0sq = nvis::inner(__v0, __v0);
    double u = v0sq / (v0sq - nvis::inner(__v0, __v1));
    nvis::vec2 x = (1 - u) * x0 + u * x1;
    return rotation_angle_linear_predictor(x0, x, __v0, __v0, true, false,
                                           map, period, dtheta, dx, depth + 1,
                                           steps, conservative) +
           rotation_angle_linear_predictor(x, x1, __v1, __v1, false, true,
                                           map, period, dtheta, dx, depth + 1,
                                           steps, conservative);
}

template<typename RHS, typename J>
double rotation_angle_cubic_predictor(const nvis::vec2& x0, const nvis::vec2& x1,
                                      // v0 and v1 are the values associated with x0 and x1...
                                      const nvis::vec2& v0, const nvis::vec2& v1,
                                      // J0 and J1 are the associated Jacobian matrices
                                      const nvis::mat2& J0, const nvis::mat2& J1,
                                      // ...provided their respective flag valid? is true
                                      bool valid0, bool valid1,
                                      const RHS& rhs, const J& jacobian,
                                      double dtheta, double dx,
                                      const metric_type& metric,
                                      unsigned int depth,
                                      std::list<step_type>& steps,
                                      bool conservative = false) const
{
    // needed for internal consistency
    const double large_angle = 0.75 * M_PI;

    nvis::vec2 __v0, __v1;
    __v0 = (valid0 ? v0 : rhs(x0));
    __v1 = (valid1 ? v1 : rhs(x1));
    nvis::mat2 __J0, __J1;
    __J0 = (valid0 ? J0 : jacobian(x0));
    __J1 = (valid1 ? J1 : jacobian(x1));

    // easy cases first
    double theta = angle(__v0, __v1);
    if (fabs(theta) < std::max(large_angle, dtheta)) return theta;
    else if (nvis::norm(x1 - x0) < dx) {
        if (conservative) {
            throw degenerate_point_exception(0.5*(x0 + x1));
        }
        return theta;
    }

    // use cubic hermite model to pick next sample
    // v(x) = v0 + j0*x + (3*(v1-v0)-2*j0-j1)*x^2 + (j1-j0-2*(v1-v0))*x^3
    nvis::vec2 e = x1 - x0;
    nvis::vec2 j0 = __J0 * e;
    nvis::vec2 j1 = __J1 * e;
    // solve for v(x).v0 = 0
    nvis::vec2 a0 = nvis::inner(__v0, __v0);
    nvis::vec2 a1 = nvis::inner(j0, __v0);
    nvis::vec2 a2 = nvis::inner(3.*(__v1 - __v0) - 2.*j0 - j1, __v0);
    nvis::vec2 a3 = nvis::inner(j1 - j0 - 2.*(__v1 - __v0), __v0);

    std::complex<double> u[3];
    int nbroots = xavier::cubic_equation(a3, a2, a1, a0, u);
    for (int i = 0 ; i < nbroots ; ++i) {
        if (u[i].imag()) continue;
        else if (u[i].real() < 0 || u[1].real() > 1) continue;

        nvis::vec2 x = (1 - u) * x0 + u * x1;
        return rotation_angle_cubic_predictor(x0, x, __v0, __v0, __J0, __J0, true, false,
                                              rhs, jacobian, dtheta, dx, depth + 1,
                                              steps, conservative) +
               rotation_angle_cubic_predictor(x, x1, __v1, __v1, __J1, __J1, false, true,
                                              rhs, jacobian, dtheta, dx, depth + 1,
                                              steps, conservative);
    }
}
}


















