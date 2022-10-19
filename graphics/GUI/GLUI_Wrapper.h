#ifndef GLUI_WRAPPER_H
#define GLUI_WRAPPER_H

#ifndef GLUT_NO_LIB_PRAGMA
#define GLUT_NO_LIB_PRAGMA
#endif
#ifndef GLUI_NO_LIB_PRAGMA
#define GLUI_NO_LIB_PRAGMA
#endif

#include <glut.h>
#include <glui.h>
#include <vector>
#include <string>
#include <math.h>

#if defined(_WIN64)  // note: VS syntax highlighting may grey out below even though it is defined and will be used
#pragma comment(lib, "glut64.lib")

#if defined(_DEBUG)
#pragma comment(lib, "glui64D.lib")
#else
#pragma comment(lib, "glui64.lib")
#endif
#elif defined(_WIN32)
#pragma comment(lib, "glut32.lib")

#if defined(_DEBUG)
#pragma comment(lib, "glui32D.lib")
#else
#pragma comment(lib, "glui32.lib")
#endif
#endif

using namespace std;

// to do: be able to use unsigned int, size_t types with spinners

namespace GLUI_Wrapper
{
class GLUI_Widget_Wrapper
{
public:
	void(*callbackPtr)(int);
	void(*wrapperCallbackPtr)(int);
	int callbackValue;

	inline void doExternalCallback() {
		if (callbackPtr != NULL && callbackPtr != wrapperCallbackPtr) callbackPtr(callbackValue);
	}

	inline void setCallbacks(int wrapperCallbackValue, void(*wrapperCallback)(int), int externalCallbackValue, void(*externalCallback)(int)) {
		callbackPtr = externalCallback;
		callbackValue = externalCallbackValue;
		wrapperCallbackPtr = wrapperCallback;
	}
};

// -------------------------------------------------------------------------- //
// --------------------------------- Panels --------------------------------- //
// -------------------------------------------------------------------------- //
static GLUI_Panel* addPanel(GLUI *glui, GLUI_Panel *panel, string name, bool showOutline = true)
{
	GLUI_Panel *newPanel = NULL;
	if (panel == NULL) newPanel = glui->add_panel((char*)name.c_str(), showOutline);
	else newPanel = glui->add_panel_to_panel(panel, (char*)name.c_str(), showOutline);
	return newPanel;
}

static GLUI_Rollout* addRollout(GLUI *glui, GLUI_Panel *panel, string name, bool open = false)
{
	GLUI_Rollout *newPanel = NULL;
	if (panel == NULL) newPanel = glui->add_rollout((char*)name.c_str(), open);
	else newPanel = glui->add_rollout_to_panel(panel, (char*)name.c_str(), open);
	return newPanel;
}

static void removePanel(GLUI_Panel *panel)
{
	if (panel != NULL) {
		panel->unlink();
		panel->disable();
		panel = NULL;
	}
}

// ------------------------------------------------------------------------------ //
// --------------------------------- Checkboxes --------------------------------- //
// ------------------------------------------------------------------------------ //

class GLUI_Checkbox_Wrapper : public GLUI_Widget_Wrapper
{
public:
	GLUI_Checkbox_Wrapper(GLUI *glui, GLUI_Panel *panel, string name, bool *value, int wrapperCallbackValue, void(*wrapperCallback)(int),
	                      int externalCallbackValue = -1, void(*externalCallback)(int) = NULL) {

		boolPtr = value;

		int *intPtr = new int(0);
		if (*boolPtr) *intPtr = 1;
		if (panel == NULL) checkbox = glui->add_checkbox((char*)name.c_str(), intPtr, wrapperCallbackValue, wrapperCallback);
		else checkbox = glui->add_checkbox_to_panel(panel, (char*)name.c_str(), intPtr, wrapperCallbackValue, wrapperCallback);

		setCallbacks(wrapperCallbackValue, wrapperCallback, externalCallbackValue, externalCallback);
	}

	GLUI_Checkbox *checkbox;
	bool *boolPtr;
};

static std::vector<GLUI_Checkbox_Wrapper> checkboxes;

static void checkboxWrapperCallback(int control)
{
	GLUI_Checkbox_Wrapper *checkbox = &checkboxes[control];
	int value = checkbox->checkbox->get_int_val();
	bool boolVal = true;
	if (value == 0) boolVal = false;
	*checkbox->boolPtr = boolVal;

	//for (unsigned int i=0; i<checkboxes.size(); i++) {
	//	if (checkboxes[i].boolPtr == checkbox->boolPtr) *checkboxes[i].boolPtr = boolVal;
	//}

	//if (*checkbox->boolPtr) printf("%s true\n", string(checkbox->checkbox->name).c_str());
	//else printf("%s false\n", string(checkbox->checkbox->name).c_str());

	checkbox->doExternalCallback();
}

static GLUI_Checkbox* addCheckbox(GLUI *glui, GLUI_Panel *panel, string name, bool *value, int callbackValue = -1, void(*callback)(int) = NULL)
{
	checkboxes.push_back(GLUI_Checkbox_Wrapper(glui, panel, name, value, (int)checkboxes.size(), checkboxWrapperCallback, callbackValue, callback));
	return checkboxes[checkboxes.size()-1].checkbox;
}

// ---------------------------------------------------------------------------- //
// --------------------------------- Spinners --------------------------------- //
// ---------------------------------------------------------------------------- //

class GLUI_Spinner_Wrapper : public GLUI_Widget_Wrapper
{
public:
	GLUI_Spinner_Wrapper(GLUI *glui, GLUI_Panel *panel, string name, double *value, int wrapperCallbackValue, void(*wrapperCallback)(int),
	                     int externalCallbackValue = -1, void(*externalCallback)(int) = NULL) {

		setValue(NULL, NULL, NULL, value);

		float *ptr = new float(*doublePtr);
		if (panel == NULL) spinner = glui->add_spinner((char*)name.c_str(), GLUI_SPINNER_FLOAT, ptr, wrapperCallbackValue, wrapperCallback);
		else spinner = glui->add_spinner_to_panel(panel, (char*)name.c_str(), GLUI_SPINNER_FLOAT, ptr, wrapperCallbackValue, wrapperCallback);

		setCallbacks(wrapperCallbackValue, wrapperCallback, externalCallbackValue, externalCallback);
	}

	GLUI_Spinner_Wrapper(GLUI *glui, GLUI_Panel *panel, string name, int *value, int wrapperCallbackValue, void(*wrapperCallback)(int),
	                     int externalCallbackValue = -1, void(*externalCallback)(int) = NULL) {

		setValue(value, NULL, NULL, NULL);

		if (panel == NULL) spinner = glui->add_spinner((char*)name.c_str(), GLUI_SPINNER_INT, intPtr, wrapperCallbackValue, wrapperCallback);
		else spinner = glui->add_spinner_to_panel(panel, (char*)name.c_str(), GLUI_SPINNER_INT, intPtr, wrapperCallbackValue, wrapperCallback);

		setCallbacks(wrapperCallbackValue, wrapperCallback, externalCallbackValue, externalCallback);
	}

	GLUI_Spinner_Wrapper(GLUI *glui, GLUI_Panel *panel, string name, float *value, int wrapperCallbackValue, void(*wrapperCallback)(int),
	                     int externalCallbackValue = -1, void(*externalCallback)(int) = NULL) {

		setValue(NULL, NULL, value, NULL);

		if (fabs(*floatPtr) < 1e-7) *floatPtr = 0;

		if (panel == NULL) spinner = glui->add_spinner((char*)name.c_str(), GLUI_SPINNER_FLOAT, floatPtr, wrapperCallbackValue, wrapperCallback);
		else spinner = glui->add_spinner_to_panel(panel, (char*)name.c_str(), GLUI_SPINNER_FLOAT, floatPtr, wrapperCallbackValue, wrapperCallback);

		setCallbacks(wrapperCallbackValue, wrapperCallback, externalCallbackValue, externalCallback);
	}

	GLUI_Spinner_Wrapper(GLUI *glui, GLUI_Panel *panel, string name, unsigned int *value, int wrapperCallbackValue, void(*wrapperCallback)(int),
	                     int externalCallbackValue = -1, void(*externalCallback)(int) = NULL) {

		setValue(NULL, value, NULL, NULL);

		int *ptr = new int(*uintPtr);
		if (panel == NULL) spinner = glui->add_spinner((char*)name.c_str(), GLUI_SPINNER_INT, ptr, wrapperCallbackValue, wrapperCallback);
		else spinner = glui->add_spinner_to_panel(panel, (char*)name.c_str(), GLUI_SPINNER_INT, ptr, wrapperCallbackValue, wrapperCallback);

		setCallbacks(wrapperCallbackValue, wrapperCallback, externalCallbackValue, externalCallback);
	}

	void setValue(int *intVal, unsigned int *uintVal, float *floatVal, double *doubleVal) {
		intPtr = intVal;
		uintPtr = uintVal;
		floatPtr = floatVal;
		doublePtr = doubleVal;
	}

	inline void clamp(float low, float high) {
		if (low < high) {
			if (intPtr != NULL || uintPtr != NULL) {
				if (uintPtr != NULL) if (low < 0) low = 0;
				spinner->set_int_limits(low, high);
			}
			else spinner->set_float_limits(low, high);
		}
	}

	GLUI_Spinner *spinner;
	double *doublePtr;
	float *floatPtr;
	int *intPtr;
	unsigned int *uintPtr;
};

static std::vector<GLUI_Spinner_Wrapper> spinners;

static void spinnerWrapperCallback(int control)
{
	GLUI_Spinner_Wrapper *spinner = &spinners[control];
	if (spinner->doublePtr != NULL) *spinner->doublePtr = spinner->spinner->get_float_val();
	if (spinner->uintPtr != NULL) *spinner->uintPtr = spinner->spinner->get_int_val();

	//if (spinner->doublePtr != NULL) printf("%s %f\n", string(spinner->spinner->name).c_str(), *spinner->doublePtr);
	//if (spinner->floatPtr != NULL) printf("%s %f\n", string(spinner->spinner->name).c_str(), *spinner->floatPtr);
	//if (spinner->intPtr != NULL) printf("%s %d\n", string(spinner->spinner->name).c_str(), *spinner->intPtr);
	//if (spinner->uintPtr != NULL) printf("%s %d\n", string(spinner->spinner->name).c_str(), *spinner->uintPtr);

	spinner->doExternalCallback();
}

template<class T>
static GLUI_Spinner* addSpinner(GLUI *glui, GLUI_Panel *panel, string name, T *value, int callbackValue = -1, void(*callback)(int) = NULL)
{
	spinners.push_back(GLUI_Spinner_Wrapper(glui, panel, name, value, (int)spinners.size(), spinnerWrapperCallback, callbackValue, callback));
	return spinners[spinners.size()-1].spinner;
}

enum { VECTOR_SPINNER_REGULAR = 0, VECTOR_SPINNER_VERTICAL, VECTOR_SPINNER_HORIZONTAL };

#if 0
struct SpinnerArray {
	vector<GLUI_Spinner*> spinners;
};

template <class T>
static SpinnerArray addSpinner(GLUI *glui, GLUI_Panel *panel, string name, T *value, unsigned int arraySize, string extension = "", int type = VECTOR_SPINNER_REGULAR, int callbackValue = -1, void(*callback)(int) = NULL)
{
	vector<string> extensions = tokenfy(extension);
	if (extensions.size() > 0 && extensions.size() != arraySize) {
		printf("Warning: Extension string for vector spinner %s does not match array size, not adding them\n", name.c_str());
		extensions.clear();
	}

	GLUI_Panel *vecPanel = panel;
	if (type == VECTOR_SPINNER_VERTICAL || type == VECTOR_SPINNER_HORIZONTAL) vecPanel = addPanel(glui, panel, name, 1);

	SpinnerArray spinnerArray;
	for (unsigned int i = 0; i < arraySize; i++) {
		string fullName = name;
		if (extensions.size() != 0) {
			if (type == VECTOR_SPINNER_REGULAR) fullName += "." + extensions[i];
			else fullName = extensions[i];
		}

		spinnerArray.spinners.push_back(addSpinner(glui, vecPanel, fullName, &value[i], callbackValue, callback));
		if (i != arraySize - 1 && type == VECTOR_SPINNER_HORIZONTAL) glui->add_column_to_panel(vecPanel, 0);
	}
	return spinnerArray;
}
#endif

static void clampSpinner(float lowVal, float highVal, GLUI_Spinner *spinner)
{
	for (unsigned int i = 0; i < spinners.size(); i++) if (spinners[i].spinner == spinner) {
			spinners[i].clamp(lowVal, highVal);
			break;
		}
}
template <class T> static void clampLastAddedSpinners(T *lowVals, T *highVals, unsigned int numberOfSpinners)
{
	for (unsigned int i = 0; i < numberOfSpinners; i++) clampSpinner(lowVals[numberOfSpinners-1-i], highVals[numberOfSpinners-1-i], spinners[spinners.size()-1 - i].spinner);
}
static void clampLastAddedSpinners(float lowVal, float highVal, unsigned int numberOfSpinners)
{
	for (unsigned int i = 0; i < numberOfSpinners; i++) clampSpinner(lowVal, highVal, spinners[spinners.size()-1 - i].spinner);
}
static void clampLastAddedSpinner(float lowVal, float highVal)
{
	clampLastAddedSpinners(lowVal, highVal, 1);
}

// ------------------------------------------------------------------------------ //
// --------------------------------- List Boxes --------------------------------- //
// ------------------------------------------------------------------------------ //

class GLUI_Listbox_Wrapper : public GLUI_Widget_Wrapper
{
public:
	GLUI_Listbox_Wrapper(GLUI *glui, GLUI_Panel *panel, string name, int *value, int wrapperCallbackValue, void(*wrapperCallback)(int),
	                     int externalCallbackValue = -1, void(*externalCallback)(int) = NULL) {

		intPtr = value;
		if (intPtr == NULL) intPtr = new int(0);

		if (panel == NULL) listbox = glui->add_listbox((char*)name.c_str(), intPtr, wrapperCallbackValue, wrapperCallback);
		else listbox = glui->add_listbox_to_panel(panel, (char*)name.c_str(), intPtr, wrapperCallbackValue, wrapperCallback);

		setCallbacks(wrapperCallbackValue, wrapperCallback, externalCallbackValue, externalCallback);
	}

	GLUI_Listbox *listbox;
	int *intPtr;
};

static std::vector<GLUI_Listbox_Wrapper> listboxes;

static void listboxWrapperCallback(int control)
{
	GLUI_Listbox_Wrapper *listbox = &listboxes[control];

	listbox->doExternalCallback();
}

static GLUI_Listbox* addListbox(GLUI *glui, GLUI_Panel *panel, string name, int *val = NULL, int callbackValue = -1, void(*callback)(int) = NULL)
{
	listboxes.push_back(GLUI_Listbox_Wrapper(glui, panel, name, val, (int)listboxes.size(), listboxWrapperCallback, callbackValue, callback));
	return listboxes[listboxes.size()-1].listbox;
}

static void addItemToListbox(GLUI_Listbox *listbox, string text, int id)
{
	listbox->add_item(id, (char*)text.c_str());
}
static void addItemToLastListbox(string text, int id)
{
	addItemToListbox(listboxes[listboxes.size()-1].listbox, text, id);
}

static void setListboxValue(GLUI_Listbox *listbox, int val)
{
	listbox->set_int_val(val);
}
static void setLastListboxValue(int val)
{
	setListboxValue(listboxes[listboxes.size()-1].listbox, val);
}

// ------------------------------------------------------------------------------ //
// -------------------------------- Radio Buttons ------------------------------- //
// ------------------------------------------------------------------------------ //

class GLUI_RadioButtonGroup_Wrapper : public GLUI_Widget_Wrapper
{
public:
	GLUI_RadioButtonGroup_Wrapper(GLUI *glui, GLUI_Panel *panel, int *value, int wrapperCallbackValue, void(*wrapperCallback)(int),
	                              int externalCallbackValue = -1, void(*externalCallback)(int) = NULL) {

		intPtr = value;
		if (intPtr == NULL) intPtr = new int(0);

		if (panel == NULL) radioGroup = glui->add_radiogroup(intPtr, wrapperCallbackValue, wrapperCallback);
		else radioGroup = glui->add_radiogroup_to_panel(panel, intPtr, wrapperCallbackValue, wrapperCallback);

		setCallbacks(wrapperCallbackValue, wrapperCallback, externalCallbackValue, externalCallback);
	}

	GLUI_RadioGroup *radioGroup;
	int *intPtr;

	std::vector<int> buttonIDs;
};

static std::vector<GLUI_RadioButtonGroup_Wrapper> radioButtonGroups;

static void radioButtonGroupWrapperCallback(int control)
{
	GLUI_RadioButtonGroup_Wrapper *radioGroup = &radioButtonGroups[control];

	if (radioGroup->radioGroup->get_int_val() >= (int)radioGroup->buttonIDs.size()) *radioGroup->intPtr = 0;
	else *radioGroup->intPtr = radioGroup->buttonIDs[radioGroup->radioGroup->get_int_val()];

	radioGroup->doExternalCallback();
}

static GLUI_RadioGroup* addRadioGroup(GLUI *glui, GLUI_Panel *panel, int *val = NULL, int callbackValue = -1, void(*callback)(int) = NULL)
{
	radioButtonGroups.push_back(GLUI_RadioButtonGroup_Wrapper(glui, panel, val, (int)radioButtonGroups.size(), radioButtonGroupWrapperCallback, callbackValue, callback));
	return radioButtonGroups[radioButtonGroups.size()-1].radioGroup;
}

static void addButtonToLastRadioGroup(GLUI *glui, string text, int id)
{
	GLUI_RadioButtonGroup_Wrapper *group = &radioButtonGroups[radioButtonGroups.size()-1];
	glui->add_radiobutton_to_group(group->radioGroup, text.c_str());
	group->buttonIDs.push_back(id);
}

static void updateLastRadioGroup()
{
	GLUI_RadioButtonGroup_Wrapper *radioGroup = &radioButtonGroups[radioButtonGroups.size()-1];

	for (unsigned int i = 0; i < radioGroup->buttonIDs.size(); i++) {
		if (radioGroup->buttonIDs[i] == *radioGroup->intPtr) {
			radioGroup->radioGroup->set_int_val(i);
			break;
		}
	}
}

// ------------------------------------------------------------------------------ //
// --------------------------------- Scrollbars --------------------------------- //
// ------------------------------------------------------------------------------ //

// since C++ only allows integer values used with #if, you need to edit glui.h so that the original GLUI_VERSION is multiplied by 100
// #if GLUI_VERSION >= 230
class GLUI_Scrollbar_Wrapper : public GLUI_Widget_Wrapper
{
public:
	GLUI_Scrollbar_Wrapper(GLUI *glui, GLUI_Panel *panel, string name, double *value, int orientation, int wrapperCallbackValue, void(*wrapperCallback)(int),
	                       int externalCallbackValue = -1, void(*externalCallback)(int) = NULL) {

		setValue(NULL, NULL, value);
		float *ptr = new float(*doublePtr);

		GLUI_StaticText *text = glui->add_statictext_to_panel(panel, (char*)name.c_str());
		text->set_w(1);	// make sure the text is as close as possible to the slider
		glui->add_column_to_panel(panel, false);
		scrollbar = new GLUI_Scrollbar(panel, (char*)name.c_str(), orientation, ptr, wrapperCallbackValue, wrapperCallback);

		setCallbacks(wrapperCallbackValue, wrapperCallback, externalCallbackValue, externalCallback);
	}

	GLUI_Scrollbar_Wrapper(GLUI *glui, GLUI_Panel *panel, string name, int *value, int orientation, int wrapperCallbackValue, void(*wrapperCallback)(int),
	                       int externalCallbackValue = -1, void(*externalCallback)(int) = NULL) {

		setValue(value, NULL, NULL);

		// we want to create a panel so that the text is grouped nicely with the bar
		if (panel == NULL) panel = glui->add_panel("", 0);
		else panel = glui->add_panel_to_panel(panel, "", 0);

		GLUI_StaticText *text = glui->add_statictext_to_panel(panel, (char*)name.c_str());
		text->set_w(1);	// make sure the text is as close as possible to the slider
		glui->add_column_to_panel(panel, false);
		scrollbar = new GLUI_Scrollbar(panel, (char*)name.c_str(), orientation, intPtr, wrapperCallbackValue, wrapperCallback);

		setCallbacks(wrapperCallbackValue, wrapperCallback, externalCallbackValue, externalCallback);
	}

	GLUI_Scrollbar_Wrapper(GLUI *glui, GLUI_Panel *panel, string name, float *value, int orientation, int wrapperCallbackValue, void(*wrapperCallback)(int),
	                       int externalCallbackValue = -1, void(*externalCallback)(int) = NULL) {

		setValue(NULL, value, NULL);

		if (fabs(*floatPtr) < 1e-7) *floatPtr = 0;

		// we want to create a panel so that the text is grouped nicely with the bar
		if (panel == NULL) panel = glui->add_panel("", 0);
		else panel = glui->add_panel_to_panel(panel, "", 0);

		GLUI_StaticText *text = glui->add_statictext_to_panel(panel, (char*)name.c_str());
		text->set_w(1);	// make sure the text is as close as possible to the slider
		glui->add_column_to_panel(panel, false);
		scrollbar = new GLUI_Scrollbar(panel, (char*)name.c_str(), orientation, floatPtr, wrapperCallbackValue, wrapperCallback);

		setCallbacks(wrapperCallbackValue, wrapperCallback, externalCallbackValue, externalCallback);
	}

	void setValue(int *intVal, float *floatVal, double *doubleVal) {
		intPtr = intVal;
		floatPtr = floatVal;
		doublePtr = doubleVal;
	}

	inline void clamp(float low, float high) {
		if (low < high) {
			if (intPtr != NULL) scrollbar->set_int_limits(low, high);
			else scrollbar->set_float_limits(low, high);
		}
	}

	GLUI_Scrollbar *scrollbar;
	double *doublePtr;
	float *floatPtr;
	int *intPtr;
};

static std::vector<GLUI_Scrollbar_Wrapper> scrollbars;

static void scrollbarWrapperCallback(int control)
{
	GLUI_Scrollbar_Wrapper *scrollbar = &scrollbars[control];
	if (scrollbar->doublePtr != NULL) *scrollbar->doublePtr = scrollbar->scrollbar->get_float_val();
	scrollbar->doExternalCallback();
}

//orientation = GLUI_SCROLL_HORIZONTAL or GLUI_SCROLL_VERTICAL
template<class T>
static GLUI_Scrollbar* addScrollbar(GLUI *glui, GLUI_Panel *panel, string name, T *value, int orientation, int callbackValue = -1, void(*callback)(int) = NULL)
{
	scrollbars.push_back(GLUI_Scrollbar_Wrapper(glui, panel, name, value, orientation, scrollbars.size(), scrollbarWrapperCallback, callbackValue, callback));
	return scrollbars[scrollbars.size()-1].scrollbar;
}

struct ScrollbarArray {
	std::vector<GLUI_Scrollbar*> scrollbars;
};

#if 0
template <class T>
static ScrollbarArray addScrollbar(GLUI *glui, GLUI_Panel *panel, string name, T *value, unsigned int arraySize, int orientation, string extension = "", int callbackValue = -1, void(*callback)(int) = NULL)
{
	vector<string> extensions = tokenfy(extension);
	if (extensions.size() > 0 && extensions.size() != arraySize) {
		printf("Warning: Extension string for vector spinner %s does not match array size, not adding them\n", name.c_str());
		extensions.clear();
	}

	GLUI_Panel *vecPanel = addPanel(glui, panel, name, 1);

	ScrollbarArray scrollbarArray;
	for (unsigned int i = 0; i < arraySize; i++) {
		string fullName = "";
		if (extensions.size() != 0) fullName = extensions[i];

		scrollbarArray.scrollbars.push_back(addScrollbar(glui, vecPanel, fullName, &value[i], orientation, callbackValue, callback));
		if (i != arraySize - 1) glui->add_column_to_panel(vecPanel, 0);
	}
	return scrollbarArray;
}
#endif

static void clampScrollbar(float lowVal, float highVal, GLUI_Scrollbar *scrollbar)
{
	for (unsigned int i = 0; i < scrollbars.size(); i++) if (scrollbars[i].scrollbar == scrollbar) {
			scrollbars[i].clamp(lowVal, highVal);
			break;
		}
}
static void clampLastAddedScrollbars(float lowVal, float highVal, unsigned int numberOfScrollbars)
{
	for (unsigned int i = 0; i < numberOfScrollbars; i++) clampScrollbar(lowVal, highVal, scrollbars[scrollbars.size()-1 - i].scrollbar);
}
static void clampLastAddedScrollbar(float lowVal, float highVal)
{
	clampLastAddedScrollbars(lowVal, highVal, 1);
}

// #endif

// ------------------------------------------------------------------------------ //
// --------------------------------- Buttons ------------------------------------ //
// ------------------------------------------------------------------------------ //

class GLUI_Button_Wrapper : public GLUI_Widget_Wrapper
{
public:
	GLUI_Button_Wrapper(GLUI *glui, GLUI_Panel *panel, string positiveName, string negativeName, int wrapperCallbackValue, void(*wrapperCallback)(int),
	                    bool *value = NULL, int externalCallbackValue = -1, void(*externalCallback)(int) = NULL) {

		if (value != NULL) boolPtr = value;
		else boolPtr = new bool(true);

		posName = positiveName;
		negName = negativeName;

		string name = posName;
		if (boolPtr != NULL) if (!*boolPtr) name = negName;

		if (panel == NULL) button = glui->add_button((char*)name.c_str(), wrapperCallbackValue, wrapperCallback);
		else button = glui->add_button_to_panel(panel, (char*)name.c_str(), wrapperCallbackValue, wrapperCallback);

		setCallbacks(wrapperCallbackValue, wrapperCallback, externalCallbackValue, externalCallback);
	}

	GLUI_Button *button;
	bool *boolPtr;
	string posName, negName;
};

static std::vector<GLUI_Button_Wrapper> buttons;

static void buttonWrapperCallback(int control)
{
	GLUI_Button_Wrapper *button = &buttons[control];

	*button->boolPtr = !*button->boolPtr;
	if (*button->boolPtr) button->button->set_name((char*)button->posName.c_str());
	else button->button->set_name((char*)button->negName.c_str());

	button->doExternalCallback();
}

static GLUI_Button* addButton(GLUI *glui, GLUI_Panel *panel, string name, int callbackValue = -1, void(*callback)(int) = NULL)
{
	buttons.push_back(GLUI_Button_Wrapper(glui, panel, name, name, (int)buttons.size(), buttonWrapperCallback, NULL, callbackValue, callback));
	return buttons[buttons.size()-1].button;
}

static GLUI_Button* addButton(GLUI *glui, GLUI_Panel *panel, string name, bool *value = NULL, int callbackValue = -1, void(*callback)(int) = NULL)
{
	buttons.push_back(GLUI_Button_Wrapper(glui, panel, name, name, (int)buttons.size(), buttonWrapperCallback, value, callbackValue, callback));
	return buttons[buttons.size()-1].button;
}

static GLUI_Button* addButton(GLUI *glui, GLUI_Panel *panel, string posName, string negName, bool *value, int callbackValue = -1, void(*callback)(int) = NULL)
{
	buttons.push_back(GLUI_Button_Wrapper(glui, panel, posName, negName, (int)buttons.size(), buttonWrapperCallback, value, callbackValue, callback));
	return buttons[buttons.size()-1].button;
}

// ------------------------------------------------------------------------ //
// --------------------------------- Misc --------------------------------- //
// ------------------------------------------------------------------------ //
static void syncLive()
{
	for (unsigned int i = 0; i < checkboxes.size(); i++) {
		if (*checkboxes[i].boolPtr) checkboxes[i].checkbox->set_int_val(1);
		else checkboxes[i].checkbox->set_int_val(0);
	}

	for (unsigned int i = 0; i < spinners.size(); i++) {
		if (spinners[i].doublePtr != NULL) spinners[i].spinner->set_float_val((float)*spinners[i].doublePtr);
		if (spinners[i].floatPtr != NULL) spinners[i].spinner->set_float_val(*spinners[i].floatPtr);
		if (spinners[i].intPtr != NULL) spinners[i].spinner->set_int_val(*spinners[i].intPtr);
		if (spinners[i].uintPtr != NULL) spinners[i].spinner->set_int_val(*spinners[i].uintPtr);
	}

	for (unsigned int i = 0; i < listboxes.size(); i++) {
		listboxes[i].listbox->set_int_val(*listboxes[i].intPtr);
	}

// #if GLUI_VERSION >= 230

	for (unsigned int i = 0; i < scrollbars.size(); i++) {
		if (scrollbars[i].doublePtr != NULL) scrollbars[i].scrollbar->set_float_val((float)*scrollbars[i].doublePtr);
		if (scrollbars[i].floatPtr != NULL) scrollbars[i].scrollbar->set_float_val(*scrollbars[i].floatPtr);
		if (scrollbars[i].intPtr != NULL) scrollbars[i].scrollbar->set_int_val(*scrollbars[i].intPtr);
	}

// #endif
}
}

#endif
