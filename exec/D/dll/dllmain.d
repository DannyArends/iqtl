/**********************************************************************
 * src/D/dll/dllmain.d
 *
 * copyright (c) 2012 Danny Arends
 * last modified Jan, 2012
 * first written Jan, 2012
 **********************************************************************/

version(windows){

  import std.c.windows.windows;
  import core.sys.windows.dll;

  __gshared HINSTANCE g_hInst;

  extern(Windows)
  BOOL DllMain(HINSTANCE hInstance, ULONG ulReason, LPVOID pvReserved){
    switch (ulReason){
      case DLL_PROCESS_ATTACH:
        g_hInst = hInstance;
        dll_process_attach( hInstance, true );
      break;
      case DLL_PROCESS_DETACH:
        dll_process_detach( hInstance, true );
      break;
      case DLL_THREAD_ATTACH:
        dll_thread_attach( true, true );
      break;
      case DLL_THREAD_DETACH:
        dll_thread_detach( true, true );
      break;
      default:
      break;
    }
    return true;
  }

}