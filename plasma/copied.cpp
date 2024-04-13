#include <windows.h>
#include <math.h>

char progname[]="Cute plasma";
char SINTAB[256];

// forward declaration:
LRESULT CALLBACK winproc(HWND hWnd, UINT uMsg, WPARAM wParam, LPARAM lParam);

int WINAPI WinMain(HINSTANCE hInst, HINSTANCE hPrevInst, LPSTR lpszCmdLine, int nCmdShow)
{
  WNDCLASSEX winclass;
  HWND hWnd;
  MSG msg;
  int i;
  for (i=0;i<256;i++)
    SINTAB[i]=sin(((i+1)*3.14159265359)/128)*127+128;

  winclass.cbSize=sizeof(WNDCLASSEX);
  winclass.style=CS_DBLCLKS;
  winclass.lpfnWndProc=&winproc;
  winclass.cbClsExtra=0;
  winclass.cbWndExtra=0;
  winclass.hInstance=hInst;
  winclass.hIcon=LoadIcon(NULL,IDI_WINLOGO);
  winclass.hCursor=LoadCursor(NULL,IDC_NO);
  winclass.hbrBackground=NULL;
  winclass.lpszMenuName=NULL;
  winclass.lpszClassName=progname;
  winclass.hIconSm=NULL;

  if (!RegisterClassEx(&winclass))
    return 0;
  hWnd=CreateWindow(
    progname,
    progname,
    WS_SYSMENU|WS_CAPTION|WS_BORDER|WS_OVERLAPPED|WS_VISIBLE|WS_MINIMIZEBOX,
    CW_USEDEFAULT,
    0,
    320+2,
    200+16+2,
    NULL,
    NULL,
    hInst,
    NULL);
  ShowWindow(hWnd,nCmdShow);
  UpdateWindow(hWnd);
  while (GetMessage(&msg,NULL,0,0))
  {
    TranslateMessage(&msg);
    DispatchMessage(&msg);
  }
  return (msg.wParam);
}

HDC pDC;
HBITMAP old;
HBITMAP ourbitmap;
int * framebuf;

void render_effect(int tick,int * framebuf)
{
  int i,j,k;
  tick/=4;
  for (k=0,i=0;i<200;i++)
    for (j=0;j<320;j++,k++)
      *(framebuf+k)=RGB(SINTAB[(i+tick)&0xff],
        SINTAB[(j-tick)&0xff],
        SINTAB[(SINTAB[tick&0xff]+(k>>6))&0xff]);
}

void render(HDC hDC)
{
  render_effect(GetTickCount(),framebuf);
  BitBlt(hDC,0,0,320,200,pDC,0,0,SRCCOPY);
}

void deinit_framebuf(void)
{
  SelectObject(pDC,old);
  DeleteDC(pDC);
  DeleteObject(ourbitmap);
}

void init_framebuf(void)
{
  HDC hDC;
  BITMAPINFO bitmapinfo;
  hDC=CreateCompatibleDC(NULL);
  bitmapinfo.bmiHeader.biSize=sizeof(BITMAPINFOHEADER);
  bitmapinfo.bmiHeader.biWidth=320;
  bitmapinfo.bmiHeader.biHeight=-200; /* top-down */
  bitmapinfo.bmiHeader.biPlanes=1;
  bitmapinfo.bmiHeader.biBitCount=32;
  bitmapinfo.bmiHeader.biCompression=BI_RGB;
  bitmapinfo.bmiHeader.biSizeImage=0;
  bitmapinfo.bmiHeader.biClrUsed=256;
  bitmapinfo.bmiHeader.biClrImportant=256;
  ourbitmap=CreateDIBSection(hDC,&bitmapinfo,DIB_RGB_COLORS,(void**)&framebuf,0,0);
  pDC=CreateCompatibleDC(NULL);
  old=(HBITMAP)SelectObject(pDC,ourbitmap);
  DeleteDC(hDC);
}

LRESULT CALLBACK winproc(HWND hWnd, UINT uMsg, WPARAM wParam, LPARAM lParam)
{
  HDC hDC;
  PAINTSTRUCT PtStr;
  switch (uMsg) {
  case WM_DESTROY:
    deinit_framebuf();
    PostQuitMessage(0);
    KillTimer (hWnd, 1);
    break;
  case WM_CREATE:
    SetTimer (hWnd, 1, 1, NULL);
    init_framebuf();
    break;
  case WM_TIMER:
    InvalidateRgn(hWnd,0,0);
    UpdateWindow (hWnd);
    break;
  case WM_PAINT:
    hDC=BeginPaint(hWnd,&PtStr);
    render(hDC);
    EndPaint(hWnd,&PtStr);
    break;
  default:
    return DefWindowProc (hWnd, uMsg, wParam, lParam);
    break;
  }
  return 0;
}