
var menutimeout  = 500;
var closetimer  = 0;
var ddmenuitem  = 0;
var menuopentimeout=150;
var opentimer=0;

// mouseover calls this
function menu_open(id)
{
  menu_cancelclosetime();
  if (ddmenuitem)
    menu_show(id);
  else
  {
    opentimer = window.setTimeout(function()
    {
      menu_show(id);
    }, menuopentimeout);
  }
}

// open menu
function menu_show(id)
{  
  // cancel close timer
  menu_cancelclosetime();

  // close old layer
  if(ddmenuitem) ddmenuitem.style.visibility = 'hidden';

  // get new layer and show it
  ddmenuitem = document.getElementById(id);
  ddmenuitem.style.visibility = 'visible';

}
// close showed layer
function menu_close()
{
  if(ddmenuitem) ddmenuitem.style.visibility = 'hidden';
        ddmenuitem=0;
}

// go close timer - mouseout calls this
function menu_closetime()
{
  if(opentimer)
  {
    window.clearTimeout(opentimer);
   opentimer = null;
  }
  closetimer = window.setTimeout(menu_close, menutimeout);
}

// cancel close timer
function menu_cancelclosetime()
{
  if(closetimer)
  {
    window.clearTimeout(closetimer);
    closetimer = null;
  }
}

// close layer when click-out
document.onclick = menu_close; 
