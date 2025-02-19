#pragma once

#include <vector>
#include <math/types.hpp>
#include <math/bounding_box.hpp>
#include <math/basic_math.hpp>

template <typename TensorField>
class EigenvectorField {
public:
    EigenvectorField(const TensorField& field)
        : tfield(field) {}
        
    const spurt::bbox3 bounds() const {
        return tfield.bounds();
    }
    
    spurt::vec3 evec(const spurt::vec3& x, int which) const {
        assert(which >= 0 && which < 3);
        mat3 tensor = spurt::from_dti(x.data());
        vec3 evals;
        mat3 evecs;
        spurt::sym_eigensystem(evals, evecs, tensor);
        return evecs.column(which);
    }
    
    double eval(const spurt::vec3& x, int which) const {
        assert(which >= 0 && which < 3);
        mat3 tensor = spurt::from_dti(x.data());
        vec3 evals;
        mat3 evecs;
        spurt::sym_eigensystem(evals, evecs, tensor);
        return evals[which];
    }
    
    spurt::vec4 eigen(const spurt::vec3& x, int which) const {
        assert(which >= 0 && which < 3);
        mat3 tensor = spurt::from_dti(x.data());
        vec3 evals;
        mat3 evecs;
        spurt::sym_eigensystem(evals, evecs, tensor);
        return spurt::vec4(evecs(0, which), evecs(1,which), evecs(2,which), evals[which]);
    }
    
private:
    const TensorField&  tfield;
};












