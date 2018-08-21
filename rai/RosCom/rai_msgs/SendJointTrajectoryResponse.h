// Generated by gencpp from file rai_msgs/SendJointTrajectoryResponse.msg
// DO NOT EDIT!


#ifndef RAI_MSGS_MESSAGE_SENDJOINTTRAJECTORYRESPONSE_H
#define RAI_MSGS_MESSAGE_SENDJOINTTRAJECTORYRESPONSE_H


#include <string>
#include <vector>
#include <map>

#include <ros/types.h>
#include <ros/serialization.h>
#include <ros/builtin_message_traits.h>
#include <ros/message_operations.h>


namespace rai_msgs
{
template <class ContainerAllocator>
struct SendJointTrajectoryResponse_
{
  typedef SendJointTrajectoryResponse_<ContainerAllocator> Type;

  SendJointTrajectoryResponse_()
    : success(false)  {
    }
  SendJointTrajectoryResponse_(const ContainerAllocator& _alloc)
    : success(false)  {
  (void)_alloc;
    }



   typedef uint8_t _success_type;
  _success_type success;





  typedef boost::shared_ptr< ::rai_msgs::SendJointTrajectoryResponse_<ContainerAllocator> > Ptr;
  typedef boost::shared_ptr< ::rai_msgs::SendJointTrajectoryResponse_<ContainerAllocator> const> ConstPtr;

}; // struct SendJointTrajectoryResponse_

typedef ::rai_msgs::SendJointTrajectoryResponse_<std::allocator<void> > SendJointTrajectoryResponse;

typedef boost::shared_ptr< ::rai_msgs::SendJointTrajectoryResponse > SendJointTrajectoryResponsePtr;
typedef boost::shared_ptr< ::rai_msgs::SendJointTrajectoryResponse const> SendJointTrajectoryResponseConstPtr;

// constants requiring out of line definition



template<typename ContainerAllocator>
std::ostream& operator<<(std::ostream& s, const ::rai_msgs::SendJointTrajectoryResponse_<ContainerAllocator> & v)
{
ros::message_operations::Printer< ::rai_msgs::SendJointTrajectoryResponse_<ContainerAllocator> >::stream(s, "", v);
return s;
}

} // namespace rai_msgs

namespace ros
{
namespace message_traits
{



// BOOLTRAITS {'IsFixedSize': True, 'IsMessage': True, 'HasHeader': False}
// {'geometry_msgs': ['/opt/ros/kinetic/share/geometry_msgs/msg'], 'trajectory_msgs': ['/opt/ros/kinetic/share/trajectory_msgs/msg'], 'std_msgs': ['/opt/ros/kinetic/share/std_msgs/msg'], 'rai_msgs': ['/home/mtoussai/git/rai-python/rai/rai/rai_msgs/msg']}

// !!!!!!!!!!! ['__class__', '__delattr__', '__dict__', '__doc__', '__eq__', '__format__', '__getattribute__', '__hash__', '__init__', '__module__', '__ne__', '__new__', '__reduce__', '__reduce_ex__', '__repr__', '__setattr__', '__sizeof__', '__str__', '__subclasshook__', '__weakref__', '_parsed_fields', 'constants', 'fields', 'full_name', 'has_header', 'header_present', 'names', 'package', 'parsed_fields', 'short_name', 'text', 'types']




template <class ContainerAllocator>
struct IsFixedSize< ::rai_msgs::SendJointTrajectoryResponse_<ContainerAllocator> >
  : TrueType
  { };

template <class ContainerAllocator>
struct IsFixedSize< ::rai_msgs::SendJointTrajectoryResponse_<ContainerAllocator> const>
  : TrueType
  { };

template <class ContainerAllocator>
struct IsMessage< ::rai_msgs::SendJointTrajectoryResponse_<ContainerAllocator> >
  : TrueType
  { };

template <class ContainerAllocator>
struct IsMessage< ::rai_msgs::SendJointTrajectoryResponse_<ContainerAllocator> const>
  : TrueType
  { };

template <class ContainerAllocator>
struct HasHeader< ::rai_msgs::SendJointTrajectoryResponse_<ContainerAllocator> >
  : FalseType
  { };

template <class ContainerAllocator>
struct HasHeader< ::rai_msgs::SendJointTrajectoryResponse_<ContainerAllocator> const>
  : FalseType
  { };


template<class ContainerAllocator>
struct MD5Sum< ::rai_msgs::SendJointTrajectoryResponse_<ContainerAllocator> >
{
  static const char* value()
  {
    return "358e233cde0c8a8bcfea4ce193f8fc15";
  }

  static const char* value(const ::rai_msgs::SendJointTrajectoryResponse_<ContainerAllocator>&) { return value(); }
  static const uint64_t static_value1 = 0x358e233cde0c8a8bULL;
  static const uint64_t static_value2 = 0xcfea4ce193f8fc15ULL;
};

template<class ContainerAllocator>
struct DataType< ::rai_msgs::SendJointTrajectoryResponse_<ContainerAllocator> >
{
  static const char* value()
  {
    return "rai_msgs/SendJointTrajectoryResponse";
  }

  static const char* value(const ::rai_msgs::SendJointTrajectoryResponse_<ContainerAllocator>&) { return value(); }
};

template<class ContainerAllocator>
struct Definition< ::rai_msgs::SendJointTrajectoryResponse_<ContainerAllocator> >
{
  static const char* value()
  {
    return "bool success\n\
\n\
";
  }

  static const char* value(const ::rai_msgs::SendJointTrajectoryResponse_<ContainerAllocator>&) { return value(); }
};

} // namespace message_traits
} // namespace ros

namespace ros
{
namespace serialization
{

  template<class ContainerAllocator> struct Serializer< ::rai_msgs::SendJointTrajectoryResponse_<ContainerAllocator> >
  {
    template<typename Stream, typename T> inline static void allInOne(Stream& stream, T m)
    {
      stream.next(m.success);
    }

    ROS_DECLARE_ALLINONE_SERIALIZER
  }; // struct SendJointTrajectoryResponse_

} // namespace serialization
} // namespace ros

namespace ros
{
namespace message_operations
{

template<class ContainerAllocator>
struct Printer< ::rai_msgs::SendJointTrajectoryResponse_<ContainerAllocator> >
{
  template<typename Stream> static void stream(Stream& s, const std::string& indent, const ::rai_msgs::SendJointTrajectoryResponse_<ContainerAllocator>& v)
  {
    s << indent << "success: ";
    Printer<uint8_t>::stream(s, indent + "  ", v.success);
  }
};

} // namespace message_operations
} // namespace ros

#endif // RAI_MSGS_MESSAGE_SENDJOINTTRAJECTORYRESPONSE_H
