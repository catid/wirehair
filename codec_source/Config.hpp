// This definition overrides CAT_BUILD_DLL below.  Neuters CAT_EXPORT macro so symbols are
// neither exported or imported.
#define CAT_NEUTER_EXPORT

// This definition changes the meaning of the CAT_EXPORT macro on Windows.  When defined,
// the CAT_EXPORT macro will export the associated symbol.  When undefined, it will import it.
#define CAT_BUILD_DLL

