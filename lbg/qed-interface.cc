
#ifdef MAKE_ENHANCED_LIGAND_TOOLS

#include <boost/python.hpp>

#include "qed-interface.hh"
#include "lidia-core/rdkit-interface.hh"


// try to get the python function, or set it to null.
void
lbg_info_t::setup_silicos_it_qed_default_func() {

   // 20230722-PE Let's give up with this (QED) Charles.
   // It's so crashy and I can't understand what the problem is.
   // I have spend more than a day and am giving up with it.
   //
   // I will spend time and effort with the GTK4 build (modern Python, modern RDKit)
   //
   silicos_it_qed_default_func    = NULL;
   silicos_it_qed_properties_func = NULL;
   return;


   // Are you here (again)?
   //
   // make sure that this works:
   // $ python
   // >>> import silicos_it.descriptors
   // >>> from silicos_it.descriptors import qed

   // and this:
   // $ ./coot-bin --no-graphics --python
   // >>> import silicos_it.descriptors
   // >>> from silicos_it.descriptors import qed

   // Build the name object
   PyObject *pName = PyUnicode_FromString("silicos_it.descriptors.qed"); // 20230513-PE merge: was just "silicos_it" - hmm.
   // Load the module object
   PyObject *pModule = PyImport_Import(pName);
   if (pModule == NULL) {
      // This happens when biscu-it is not installed (which is often(?) the case)
      // 
      std::cout << "Null pModule in get_qed() - biscu-it not installed? " << std::endl;
   } else { 
      pName = PyUnicode_FromString("silicos_it.descriptors.qed");
      pModule = PyImport_Import(pName);
      // pDict is a borrowed reference
      PyObject *pDict = PyModule_GetDict(pModule);
      if (! PyDict_Check(pDict)) {
         std::cout << "pDict is not a dict"<< std::endl;
      } else {
         // pFunc is also a borrowed reference
         PyObject *pFunc = PyDict_GetItemString(pDict, "default");
         if (! pFunc) {
            std::cout << "pFunc for default is NULL" << std::endl;
         } else {
            if (PyCallable_Check(pFunc)) {
               std::cout << "Yay - storing silicos_it_qed_default_func" << std::endl;
               silicos_it_qed_default_func = pFunc;
            } else {
               std::cout << "WARNING:: in setup_silicos_it_qed_default_func() default() function is not callable"
                         << std::endl;
            }
         }

	 pFunc = PyDict_GetItemString(pDict, "properties");
	 if (pFunc) {
	    if (PyCallable_Check(pFunc)) {
	       silicos_it_qed_properties_func = pFunc;
               std::cout << "Yay - storing silicos_it_qed_properties_func" << std::endl;
	    } else {
               std::cout << "WARNING:: in setup_silicos_it_qed_default_func() properties() function is not callable"
                         << std::endl;
	    }
	 } else {
            std::cout << "WARNING:: in setup_silicos_it_qed_default_func() properties() function is NULL"
                      << std::endl;
         }

         silicos_it_qed_pads = PyDict_GetItemString(pDict, "pads2"); // or pads1 for Gerebtzoff
      }
      Py_XDECREF(pModule);
   }
   Py_DECREF(pName);
}

// needs a properly initialized (say, via rdkit_mol_sanitize()) molecule.
double
get_qed(PyObject *silicos_it_qed_default_func, const RDKit::ROMol &rdkm) {

   double r = -1.0;
   if (silicos_it_qed_default_func) {
      try {
	 PyObject *arg_list = PyTuple_New(1);
	 PyObject *rdkit_mol_py;
	    
	 RDKit::ROMol *mol_copy_p = new RDKit::ROMol(rdkm);
	 boost::shared_ptr<RDKit::ROMol> xx(mol_copy_p);
	 boost::python::object obj(xx);
	 rdkit_mol_py = obj.ptr();
	 PyTuple_SetItem(arg_list, 0, rdkit_mol_py);
	 PyObject *result_default = PyEval_CallObject(silicos_it_qed_default_func, arg_list);
	 if (false)
	    std::cout << "DEBUG:: in get_qed() silicos_it_qed_default_func is "
		      << silicos_it_qed_default_func << std::endl;
	 if (! result_default) {
	    std::cout << "Null result from silicos_it_qed_default_func" << std::endl;
	 } else { 
	    if (PyFloat_Check(result_default)) { 
	       r = PyFloat_AsDouble(result_default);
	    }
	 }

	 // OK, I'm finished with mol_copy_p and obj now.  How do I
	 // get rid of them?  (Does that happen automatically?)
	    
      }
      catch (const boost::python::error_already_set &e) {
	 std::cout << "catch error_already_set exception "
		   << std::endl;
	 PyObject *type_ptr = NULL, *value_ptr = NULL, *traceback_ptr = NULL;
	 PyErr_Fetch(&type_ptr, &value_ptr, &traceback_ptr);

	 PyObject *dest = PyUnicode_FromString("object: %s\n");
	 if (type_ptr)
	    std::cout << "error: type "
		      << PyBytes_AS_STRING(PyUnicode_AsUTF8String(PyUnicode_Format(dest, type_ptr))) << std::endl;
	 if (value_ptr)
	    std::cout << "error: value "
		      << PyBytes_AS_STRING(PyUnicode_AsUTF8String(PyUnicode_Format(dest, value_ptr))) << std::endl;
	 if (traceback_ptr)
	    std::cout << "error: traceback_ptr "
		      << PyBytes_AS_STRING(PyUnicode_AsUTF8String(PyUnicode_Format(dest, traceback_ptr))) << std::endl;
	    
      } 
      catch (...) {
	 std::cout << "catch all exception" << std::endl;
      }
   }
   return r;
}

#endif

