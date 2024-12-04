#ifndef __VTK_MACROS_HPP__
#define __VTK_MACROS_HPP__

#define VTK_CREATE(type, name)                                  \
	vtkSmartPointer<type> name = vtkSmartPointer<type>::New()

#define VTK_INPUT(receiver, giver)                              \
	receiver->SetInputConnection(giver->GetOutputPort())

#define VTK_ACTOR(actor, mapper)                                \
	vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New(); \
	actor->SetMapper(mapper)

#define VTK_DRAW(mapper_type, input, basename)                         \
	vtkSmartPointer<mapper_type> basename##_mapper = vtkSmartPointer<mapper_type>::New(); \
	basename##_mapper->SetInputConnection(input->GetOutputPort());       \
	vtkSmartPointer<vtkActor> basename##_actor = vtkSmartPointer<vtkActor>::New(); \
	basename##_actor->SetMapper(basename##_mapper)

#endif