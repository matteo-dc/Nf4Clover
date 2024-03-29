#ifndef PROP_HPP
#define PROP_HPP

#ifndef EXTERN_PROP
 #define EXTERN_PROP extern
#endif

namespace jprop
{
    void set_ins();
    
    enum ins{LO};
    
    EXTERN_PROP vector<ins> ins_list;
    EXTERN_PROP int nins;
}

namespace qprop
{
    void set_ins();
    
    enum ins{LO};
    
    EXTERN_PROP vector<ins> ins_list;
    EXTERN_PROP int nins;
}

// invert the propagator
vvprop_t invert_jprop( const vvprop_t &jprop);

#undef EXTERN_PROP

#endif